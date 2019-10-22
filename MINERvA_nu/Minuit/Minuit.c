#include <iostream>
#include <fstream> // To use ifstream
#include <string>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;
using namespace TMath;

// Definition of basics constants and tables

const int mecDim = 24;
const int eventsNumber = 1e6;
const float crossSection = 1.19727e-39;
const float alpha = 0.2;

const int pT[] = {75, 150, 250, 325, 400, 475, 550, 700, 850, 1000, 1250, 1500, 2500};
const int pL[] = {1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 8000, 10000, 15000, 20000};

const int delta_pT[] = {75, 75, 100, 75, 75, 75, 75, 150, 150, 150, 250, 250, 1000}; 		// 13
const int delta_pL[] = {500, 500, 500, 500, 500, 500, 500, 1000, 2000, 2000, 5000, 5000};	// 12

double** matrixVanilla;
double** crossNuwroRaw;
double** crossDanielRaw;
double** covarianceInverse;
double** crossErrorVanilla;
double mat[156][576];
double crossError[12][13];

TMatrixD matrix(156,576);
TMatrixD crossNuwro(12,13);
TMatrixD crossDaniel(12,13);
TMatrixD crossDanielCov(12,13);
TMatrixD crossNuwroCov(12,13);
TMatrixD covariance(156,156);

vector<int> right_line{23,47,71,95,119,143,167,191,215,239,263,287,311,335,359,383,407,431,455,479,503,527,551,575};
vector<int> bottom_line{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
vector<int> diagonal_line{0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,575};
vector<int> inner{26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,51,52,53,54,55,56,57,58,59,60,61,62,63,64,
					65,66,67,68,69,70,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,101,102,103,104,105,106,107,108,
					109,110,111,112,113,114,115,116,117,118,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,
					151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,176,177,178,179,180,181,182,183,184,185,186,187,
					188,189,190,201,202,203,204,205,206,207,208,209,210,211,212,213,214,226,227,228,229,230,231,232,233,234,235,236,
					237,238,251,252,253,254,255,256,257,258,259,260,261,262,276,277,278,279,280,281,282,283,284,285,286,301,302,303,
					304,305,306,307,308,309,310,326,327,328,329,330,331,332,333,334,351,352,353,354,355,356,357,358,376,377,378,379,
					380,381,382,401,402,403,404,405,406,426,427,428,429,430,451,452,453,454,476,477,478,501,502,526};

// pT_plot = np.array([75, 150, 250, 325, 400, 475, 550, 700, 850, 1000, 1250, 1500, 2500, 3000])
// delta_pT_plot = np.array([75, 75, 100, 75, 75, 75, 75, 150, 150, 150, 250, 250, 1000, 0])

//////////////////////

double mysigmo(double x)
{
	return 100/(1+Exp(-200*x));
}

double** importArray2D(string fileName, unsigned height, unsigned width) // Importing 2D array from a file
{
  double** array2D = 0;
  array2D = new double*[height];

  ifstream inputFile(fileName);

  for (int h = 0; h < height; h++)
  {
        array2D[h] = new double[width];

        for (int w = 0; w < width; w++)
        {
              inputFile >> array2D[h][w];
        }
  }

  return array2D;
}

// penalty function added to chi2
double punishment(double *rescale)
{
	double pun = 0;

	// for bottom side of matrix
	for (int i : bottom_line){
		pun += mysigmo(-rescale[i] / Max( rescale[i + 1], Max(rescale[i + 24], rescale[i - 1])) + (1 - alpha));
    	pun += mysigmo(rescale[i]  / Min( rescale[i + 1], Min(rescale[i + 24], rescale[i - 1])) - 1/(1 - alpha));
	}
    
    // for right side of matrix
    for (int i : right_line){  
        pun += mysigmo(-rescale[i] / Max( rescale[i + 24], Max(rescale[i - 1], rescale[i - 24])) + (1 - alpha));
        pun += mysigmo(rescale[i]  / Min( rescale[i + 24], Min(rescale[i - 1], rescale[i - 24])) - 1/(1 - alpha));
    }

    // for diagonal side of matrix
    for (int i : diagonal_line){    
        pun += mysigmo(-rescale[i] / Max( rescale[i + 1], rescale[i - 24]) + (1 - alpha));
        pun += mysigmo(rescale[i]  / Min( rescale[i + 1], rescale[i - 24]) - 1/(1 - alpha));
    }

	// for the rest 
    for (int i : inner){    
        pun += mysigmo(-rescale[i] / Max( rescale[i + 1], Max(rescale[i - 1], Max(rescale[i + 24], rescale[i - 24]))) + (1 - alpha));
        pun += mysigmo(rescale[i]  / Min( rescale[i + 1], Min(rescale[i - 1], Min(rescale[i + 24], rescale[i - 24]))) - 1/(1 - alpha));
    }

    return pun;
}

// definition of chi^2 with covariance and penalty function
double chi2_cov(double *rescale){
	TMatrixD rescaling(576,1, rescale);
	TMatrixD mec(156,1);
	TMatrixD tmp(156,1);

    mec = matrix*rescaling; // rescaled MEC
        
    tmp = (crossDanielCov - crossNuwroCov - mec);

    return ((tmp*covariance)*tmp)(0,0) + punishment(rescale);  // chi^2 = tmp.cov.tmp
}

//////////////////////////////////////

int Minuit()
{
	/////////////////////////////////////////////
	// Importing matrices 

	matrixVanilla		= importArray2D("../data/matrix_LFG.dat", 					12*13, 	24*24);
	crossNuwroRaw		= importArray2D("../data/cros_total_nuwro_LFG.dat", 		12,		13);
	crossDanielRaw		= importArray2D("../data/cros_total_daniel.dat", 			12,		13);
	covarianceInverse	= importArray2D("../data/daniel_covariance_inv.dat",		12*13,	12*13);
	crossErrorVanilla	= importArray2D("../data/cros_error_daniel.dat",			12,		13);

	//////////////////////////////////////////////
	// Matrix renormalization

	for (int i = 0; i < 12; i++){
		for (int j = 0; j < 13; ++j){	
			for (int k = 0; k < 24*24; ++k){
				mat[13*i + j][k] = matrixVanilla[13*i + j][k] * crossSection * 1e6 * 12.0 / 13.0 / eventsNumber / delta_pT[j] / delta_pL[i];
			}
		}
	}

	for (int i = 0; i < 12; ++i){
		for (int j = 0; j < 13; j++){
			if (crossErrorVanilla[i][j] == 0)
				crossError[i][j] = 1;
			else
				crossError[i][j] = crossErrorVanilla[i][j];
		}
	}

	//////////////////////////////////////////////
	// Copyng arrays values to ROOT TMarixD structure

	for (int i = 0; i < 12; i++)
		for (int j = 0; j < 13; j++){
			crossNuwro(i,j) = crossNuwroRaw[i][j];
			crossDaniel(i,j) = crossDanielRaw[i][j];
		}

	for (int i = 0; i < 156; i++)
		for (int j = 0; j < 576; j++){
			matrix(i,j) = mat[i][j];
		}

	for (int i = 0; i < 156; i++)
		for (int j = 0; j < 156; j++){
			covariance(i,j) = covarianceInverse[i][j];
		}

	crossNuwroCov  = crossNuwro;
	crossDanielCov = crossDaniel;

	crossNuwroCov.T();
	crossDanielCov.T();

	crossNuwroCov.ResizeTo(1,156);
	crossDanielCov.ResizeTo(1,156);
	/////////////////////////////////////////////
		
	// crossNuwro.Transpose(crossNuwro);

	double array[576];

	fill_n(array, 576, 1);

	cout << chi2_cov(array) << endl;

	///////////////////////////////////////////////
	// Debug
	//
	// crossNuwro.Print();
	// crossDaniel.Print();
	//
	// for (int i = 0; i < 12; i++){
	// 	for (int j = 0; j < 13; j++)
	// 		cout << crossNuwroRaw[i][j] << "\t";
	// 	cout << endl; 
	// }
	////////////////////////////////////////////////


	return 0;
}
