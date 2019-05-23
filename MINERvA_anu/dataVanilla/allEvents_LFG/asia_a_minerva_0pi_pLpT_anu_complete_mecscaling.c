#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

#include "TBrowser.h"
#include "TGraph2D.h"

#include "minerva_0pi_pLpT_anu.h"

// #include "mec_scaling_april23_corr.h"
// #include "mec_scaling_minerva_anu_april13.h"

int events = 1000000;
//int events = 250;

double Mproton = 938.272;
double Mproton2 = Mproton*Mproton;

double Mneutron = 939.565;

double mpiplus=139.57;
double mpiplus2=mpiplus*mpiplus;

double mmion=105.658;
double mmion2=mmion*mmion;

double MyPi = 3.14159265358979323846;
double kosmin = cos (20.0*MyPi/180.0);

double max_nu=0;
double max_anu=0;

char napis [130];
char napis2 [130];
char napis3 [130];
char napiss [130];
char napiss2 [130];
char napiss3 [130];
char napiss4 [230];
char napiss5 [230];
char napiss6 [230];
char napisss4 [230];
char napisss5 [230];
char napisss6 [230];

const int transbins_anu = 6; 
const int longbins_anu = 10;

double xdpt_anu[transbins_anu];
double ydpt_anu[transbins_anu];
double exdpt_anu[transbins_anu];
double eydpt_anu[transbins_anu];

double ydpt_anu_norm[transbins_anu];

double trans_anu_left [transbins_anu+1] = {0, 150, 250, 400, 700, 1000, 1500};
double long_anu_left [longbins_anu+1] = {1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000, 8000, 10000, 15000};

double results_anu [longbins_anu+2] [transbins_anu+2];
double cross_results_anu [longbins_anu+1] [transbins_anu+1];
double results_anu_qel [longbins_anu+2] [transbins_anu+2];
double cross_results_anu_qel [longbins_anu+1] [transbins_anu+1];

double results_anu_res [longbins_anu+2] [transbins_anu+2];
double cross_results_anu_res [longbins_anu+1] [transbins_anu+1];

double results_anu_mec [longbins_anu+2] [transbins_anu+2];
double cross_results_anu_mec [longbins_anu+1] [transbins_anu+1];

double results_anu_dis [longbins_anu+2] [transbins_anu+2];
double cross_results_anu_dis [longbins_anu+1] [transbins_anu+1];

double cross_results_anu_norm [longbins_anu+1] [transbins_anu+1];

double q_analyser_qel [longbins_anu+2][transbins_anu+2];
double q2_analyser_qel [longbins_anu+2][transbins_anu+2];
double q_analyser_mec [longbins_anu+2][transbins_anu+2];
double q2_analyser_mec [longbins_anu+2][transbins_anu+2];
double w_analyser_qel [longbins_anu+2][transbins_anu+2];
double w2_analyser_qel [longbins_anu+2][transbins_anu+2];
double w_analyser_mec [longbins_anu+2][transbins_anu+2];
double w2_analyser_mec [longbins_anu+2][transbins_anu+2];

double q_average_qel [longbins_anu][transbins_anu];
double q_sigma_qel [longbins_anu][transbins_anu];
double w_average_qel [longbins_anu][transbins_anu];
double w_sigma_qel [longbins_anu][transbins_anu];

double q_average_mec [longbins_anu][transbins_anu];
double q_sigma_mec [longbins_anu][transbins_anu];
double w_average_mec [longbins_anu][transbins_anu];
double w_sigma_mec [longbins_anu][transbins_anu];

double norm_anu_nuwro=0;
double norm_anu_data=0;

double norm_anu_nuwro_bins [longbins_anu];
double norm_anu_data_bins [longbins_anu];

///////////////////////////////////////////////////////
// double mscaling (double momtran, double entran)
// {
//   int momnumber = int(momtran/100);
//   int ennumber = int(entran/100);
  
//   //cout<<momtran<<"  "<<momnumber<<"  "<<entran<<"  "<<ennumber<<"  "<<mec_scaling [momnumber][ennumber]<<endl;
//   //return mec_scaling [momnumber][ennumber];
//   return mec_scaling_minerva_anu [momnumber][ennumber];
//   //return 1;
// }
//////////////////////////////////////////////////////

ofstream nuwro ("1709_minerva_CC0pi_pLpT_anu_raport_mecscaling.dat");

void asia_a_daniel (TFile *Input, double table_anu [longbins_anu+2] [transbins_anu+2], 
    double table_anu_qel [longbins_anu+2] [transbins_anu+2], double table_anu_res [longbins_anu+2] [transbins_anu+2],
    double table_anu_mec [longbins_anu+2] [transbins_anu+2], double table_anu_dis [longbins_anu+2] [transbins_anu+2], 
    double table_anu_qel_q [longbins_anu+2] [transbins_anu+2], 
    double table_anu_qel_q2 [longbins_anu+2] [transbins_anu+2],
    double table_anu_qel_w [longbins_anu+2] [transbins_anu+2], 
    double table_anu_qel_w2 [longbins_anu+2] [transbins_anu+2],
    double table_anu_mec_q [longbins_anu+2] [transbins_anu+2], 
    double table_anu_mec_q2 [longbins_anu+2] [transbins_anu+2],
    double table_anu_mec_w [longbins_anu+2] [transbins_anu+2], 
    double table_anu_mec_w2 [longbins_anu+2] [transbins_anu+2] )
{
        TH2D *h2 = new TH2D ("h2", "h2", longbins_anu, long_anu_left, transbins_anu, trans_anu_left); 
	TH2D *h2qel = new TH2D ("h2qel", "h2qel", longbins_anu, long_anu_left, transbins_anu, trans_anu_left); 
	TH2D *h2res = new TH2D ("h2res", "h2res", longbins_anu, long_anu_left, transbins_anu, trans_anu_left);
	TH2D *h2mec = new TH2D ("h2mec", "h2mec", longbins_anu, long_anu_left, transbins_anu, trans_anu_left);
	TH2D *h2dis = new TH2D ("h2dis", "h2dis", longbins_anu, long_anu_left, transbins_anu, trans_anu_left);
	
	TH2D *h2qel_q = new TH2D ("h2qel_q", "h2qel_q", longbins_anu, long_anu_left, transbins_anu, trans_anu_left); 
	TH2D *h2qel_q2 = new TH2D ("h2qel_q2", "h2qel_q2", longbins_anu, long_anu_left, transbins_anu, trans_anu_left); 
	TH2D *h2qel_w = new TH2D ("h2qel_w", "h2qel_w", longbins_anu, long_anu_left, transbins_anu, trans_anu_left); 
	TH2D *h2qel_w2 = new TH2D ("h2qel_w2", "h2qel_w2", longbins_anu, long_anu_left, transbins_anu, trans_anu_left); 
	
	TH2D *h2mec_q = new TH2D ("h2mec_q", "h2mec_q", longbins_anu, long_anu_left, transbins_anu, trans_anu_left); 
	TH2D *h2mec_q2 = new TH2D ("h2mec_q2", "h2mec_q2", longbins_anu, long_anu_left, transbins_anu, trans_anu_left); 
	TH2D *h2mec_w = new TH2D ("h2mec_w", "h2mec_w", longbins_anu, long_anu_left, transbins_anu, trans_anu_left); 
	TH2D *h2mec_w2 = new TH2D ("h2mec_w2", "h2mec_w2", longbins_anu, long_anu_left, transbins_anu, trans_anu_left); 
	
	TTree *tt1 = (TTree*)Input->Get("treeout");
	event * e= new event();
	tt1->SetBranchAddress("e",&e);
	
	for( int k=0; k < events; k++ )
	{
	tt1->GetEntry(k);
	
		  if (e->fof(211, -211, 111)==0 && e->fof(321, -321, 311)==0 && e->fof(-311, 130, 310)==0)//no pions 
		  {
		    double kos = e->out[0].z/e->out[0].momentum();
		    if (kos>kosmin)
		    {
		    double muonlong = e->out[0].z;
		    double muontrans = sqrt( (e->out[0].x)*(e->out[0].x)+(e->out[0].y)*(e->out[0].y) );
		    double entrans = e->in[0].t - e->out[0].t;
		    double momtrans2 = entrans*entrans - e->q2();
		    double momtrans = sqrt(momtrans2);
		    
		    if (e->flag.mec==1)
		    h2->Fill(muonlong, muontrans/*,mscaling(momtrans,entrans)*/);
		    else
		    h2->Fill(muonlong, muontrans);
		    
		    if (e->flag.qel==1)
		    {h2qel->Fill(muonlong, muontrans);
		      h2qel_q->Fill(muonlong,muontrans,sqrt(momtrans2));
		      h2qel_q2->Fill(muonlong,muontrans,momtrans2 );
		      h2qel_w->Fill(muonlong,muontrans,entrans);
		      h2qel_w2->Fill(muonlong,muontrans,entrans*entrans );
		    }
		    if (e->flag.res==1)
		      h2res->Fill(muonlong, muontrans);
		    if (e->flag.mec==1)
		    {h2mec->Fill(muonlong, muontrans/*,mscaling(momtrans,entrans)*/);
		      h2mec_q->Fill(muonlong,muontrans,sqrt(momtrans2));
		      h2mec_q2->Fill(muonlong,muontrans,momtrans2 );
		      h2mec_w->Fill(muonlong,muontrans,entrans);
		      h2mec_w2->Fill(muonlong,muontrans,entrans*entrans );
		    }
		    if (e->flag.dis==1)
		      h2dis->Fill(muonlong, muontrans);
		    }
		  }
	}
	delete e;
	
	for (int j=0; j<longbins_anu+2; j++)
	{
	  for (int k=0; k<transbins_anu+2; k++)
	  { 
	    table_anu [j][k] = h2->GetBinContent(j,k);
	    table_anu_qel [j][k] = h2qel->GetBinContent(j,k);
	    table_anu_qel_q [j][k] = h2qel_q->GetBinContent(j,k);
	    table_anu_qel_q2 [j][k] = h2qel_q2->GetBinContent(j,k);
	    table_anu_qel_w [j][k] = h2qel_w->GetBinContent(j,k);
	    table_anu_qel_w2 [j][k] = h2qel_w2->GetBinContent(j,k);
	    table_anu_res [j][k] = h2res->GetBinContent(j,k);
	    table_anu_mec [j][k] = h2mec->GetBinContent(j,k);
	    table_anu_mec_q [j][k] = h2mec_q->GetBinContent(j,k);
	    table_anu_mec_q2 [j][k] = h2mec_q2->GetBinContent(j,k);
	    table_anu_mec_w [j][k] = h2mec_w->GetBinContent(j,k);
	    table_anu_mec_w2 [j][k] = h2mec_w2->GetBinContent(j,k);
	    table_anu_dis [j][k] = h2dis->GetBinContent(j,k);
	  }  
	}
	delete h2, h2qel, h2res, h2mec, h2dis;
}

double norm_cc (ifstream &Input)
{
string ff;
double answer=0;
double aa;
    
    while (!Input.eof())
      {
	getline (Input, ff);
	for (int kk=0; kk<10; kk++)
	{
           for (int jj=0; jj<4; jj++)
	   {Input>>aa;}
	    if (kk==0)
	   answer+=aa;
	    if (kk==2)
	   answer+=aa;
	    if (kk==4)
	   answer+=aa;
	    if (kk==6)
	   answer+=aa;
	    if (kk==8)
	   answer+=aa;
	 }
      }
return answer;
}

//////////////////////////////////////////
/////////////////////////////////////////
void asia_a_minerva_0pi_pLpT_anu_complete_mecscaling()
/////////////////////////////////////////
////////////////////////////////////////
{
ifstream Input2 ("allEvents.txt");
double cross_minerva_anu = norm_cc(Input2);

TFile *tf2 = new TFile("allEvents.root");
asia_a_daniel (tf2, results_anu, results_anu_qel, results_anu_res, results_anu_mec, results_anu_dis,
  q_analyser_qel, q2_analyser_qel, w_analyser_qel, w2_analyser_qel,
  q_analyser_mec, q2_analyser_mec, w_analyser_mec, w2_analyser_mec);
delete tf2;

nuwro<<"cross sections"<<endl;
for (int j=1; j<longbins_anu+1; j++)
	{
	  for (int k=1; k<transbins_anu+1; k++)
	  { 
	    cross_results_anu [j][k] = results_anu[j][k]/
	    events*cross_minerva_anu/(long_anu_left[j]-long_anu_left[j-1])/(trans_anu_left[k]-trans_anu_left[k-1])*1e6;
	    
	    cross_results_anu_qel [j][k] = results_anu_qel[j][k]/
	    events*cross_minerva_anu/(long_anu_left[j]-long_anu_left[j-1])/(trans_anu_left[k]-trans_anu_left[k-1])*1e6;
	    cross_results_anu_res [j][k] = results_anu_res[j][k]/
	    events*cross_minerva_anu/(long_anu_left[j]-long_anu_left[j-1])/(trans_anu_left[k]-trans_anu_left[k-1])*1e6;
	    cross_results_anu_mec [j][k] = results_anu_mec[j][k]/
	    events*cross_minerva_anu/(long_anu_left[j]-long_anu_left[j-1])/(trans_anu_left[k]-trans_anu_left[k-1])*1e6;
	    cross_results_anu_dis [j][k] = results_anu_dis[j][k]/
	    events*cross_minerva_anu/(long_anu_left[j]-long_anu_left[j-1])/(trans_anu_left[k]-trans_anu_left[k-1])*1e6;
	    nuwro<<cross_results_anu [j][k]<<"  "<<cross_results_anu_qel [j][k]<<"  "<<
	    cross_results_anu_res [j][k]<<"  "<<cross_results_anu_mec [j][k]<<"  "<<cross_results_anu_dis [j][k]<<endl;
	  }
	  nuwro<<endl;
	}
	
	for (int j=0; j<longbins_anu; j++)
	{norm_anu_nuwro_bins [j]=0;
	 norm_anu_data_bins [j]=0;
	}

////////////////////////////////////////////
//////////// q omega analyzer
///////////////////////////////////////////

nuwro<<"analysis of energy and momentum transfer"<<endl;

for (int j=1; j<longbins_anu+1; j++)
{
  for (int k=1; k<transbins_anu+1; k++)
  {
    double srednia_w_qel;
    double dyspersja_w_qel;
    double srednia_q_qel;
    double dyspersja_q_qel;
    //cout<<jj<<"  "<<k<<"  "<<results_qel[k][jj]<<"  "<<w_analyser_qel[k][jj]<<"  "<<w2_analyser_qel [k][jj]<<" next ";
    if (results_anu_qel[j][k]>0)
    {
    //nuwro<<jj<<"  "<<k<<"  "<<results_qel[k][jj]<<"  "<<w_analyser_qel[k][jj]/results_qel[k][jj]<<"  "<<w2_analyser_qel [k][jj]/results_qel[k][jj]<<" next ";
    srednia_w_qel = w_analyser_qel[j][k]/results_anu_qel[j][k];
    dyspersja_w_qel = sqrt( w2_analyser_qel [j][k]/results_anu_qel[j][k] - srednia_w_qel*srednia_w_qel );
    srednia_q_qel = q_analyser_qel[j][k]/results_anu_qel[j][k];
    dyspersja_q_qel = sqrt( q2_analyser_qel [j][k]/results_anu_qel[j][k] - srednia_q_qel*srednia_q_qel );
    }
    else
    {
      srednia_w_qel =0;
      dyspersja_w_qel=0;
      srednia_q_qel =0;
      dyspersja_q_qel=0;
    }
    w_average_qel[j-1][k-1]=srednia_w_qel;
    w_sigma_qel[j-1][k-1]=dyspersja_w_qel;
    q_average_qel[j-1][k-1]=srednia_q_qel;
    q_sigma_qel[j-1][k-1]=dyspersja_q_qel;
    nuwro<<j-1<<"  "<<k-1<<"  "<<w_average_qel[j-1][k-1]<<"  "<<w_sigma_qel[j-1][k-1]<<"  ";
    nuwro<<q_average_qel[j-1][k-1]<<"  "<<q_sigma_qel[j-1][k-1]<<endl;
  }
  nuwro<<endl;
}
cout<<endl;
nuwro<<endl;


for (int j=1; j<longbins_anu+1; j++)
{
  for (int k=1; k<transbins_anu+1; k++)
  {
    double srednia_w_mec;
    double dyspersja_w_mec;
    double srednia_q_mec;
    double dyspersja_q_mec;
    //cout<<jj<<"  "<<k<<"  "<<results_mec[k][jj]<<"  "<<w_analyser_mec[k][jj]<<"  "<<w2_analyser_mec [k][jj]<<" next ";
    if (results_anu_mec[j][k]>0)
    {
    //nuwro<<jj<<"  "<<k<<"  "<<results_mec[k][jj]<<"  "<<w_analyser_mec[k][jj]/results_mec[k][jj]<<"  "<<w2_analyser_mec [k][jj]/results_mec[k][jj]<<" next ";
    srednia_w_mec = w_analyser_mec[j][k]/results_anu_mec[j][k];
    dyspersja_w_mec = sqrt( w2_analyser_mec [j][k]/results_anu_mec[j][k] - srednia_w_mec*srednia_w_mec );
    srednia_q_mec = q_analyser_mec[j][k]/results_anu_mec[j][k];
    dyspersja_q_mec = sqrt( q2_analyser_mec [j][k]/results_anu_mec[j][k] - srednia_q_mec*srednia_q_mec );
    }
    else
    {
      srednia_w_mec =0;
      dyspersja_w_mec=0;
      srednia_q_mec =0;
      dyspersja_q_mec=0;
    }
    w_average_mec[j-1][k-1]=srednia_w_mec;
    w_sigma_mec[j-1][k-1]=dyspersja_w_mec;
    q_average_mec[j-1][k-1]=srednia_q_mec;
    q_sigma_mec[j-1][k-1]=dyspersja_q_mec;
    nuwro<<j-1<<"  "<<k-1<<"  "<<w_average_mec[j-1][k-1]<<"  "<<w_sigma_mec[j-1][k-1]<<"  ";
    nuwro<<q_average_mec[j-1][k-1]<<"  "<<q_sigma_mec[j-1][k-1]<<endl;
  }
  nuwro<<endl;
}
cout<<endl;
nuwro<<endl;

nuwro<<"end q omega analysis"<<endl;
	
///////////////////////////////
///////// chi2 ///////////////
//////////////////////////////

TMatrixD cov(60,60);
TMatrixD cov_copy(60,60);
TMatrixD cov1(6,6);
int counter1=0;
int counter2=0;
double chi2=0;
double chikwa[longbins_anu];//separate chi2 in long bins

for (int s=0; s<longbins_anu; s++)
  chikwa[s]=0;

for (int k=0; k<longbins_anu; k++)
  {
    for(int l=0; l<transbins_anu; l++)
    {
    for (int m=0; m<longbins_anu; m++)
    {
      for (int n=0; n<transbins_anu; n++)
      {//we keep binning convention from the paper
      cov [l*longbins_anu+k][n*longbins_anu+m] = minerva_anu_results_covariance [l*longbins_anu+k][n*longbins_anu+m];
     
      if ( (l*longbins_anu+k)==(n*longbins_anu+m) && cov [l*longbins_anu+k][n*longbins_anu+m] <1e-15 )
      {cov [l*longbins_anu+k][n*longbins_anu+m] = 1;//making matrix invertible
	counter1++;//control counter
      }
      cov_copy [l*longbins_anu+k][n*longbins_anu+m] = cov [l*longbins_anu+k][n*longbins_anu+m];
      }
    }
    }
  }

  cov.Invert();

for (int k=0; k<60; k++)
  {
    for (int j=0; j<60; j++)
    {
        if (k==j && cov[k][j]<1+1e-15 && cov[k][j]>1-1e-15)
	{cov[k][j]=0;
	  counter2++;//control counter
	}
    }
  }//

for (int k=0; k<longbins_anu; k++)
{
  for (int l=0; l<transbins_anu; l++)
  {
    for (int m=0; m<longbins_anu; m++)
    {
      for (int n=0; n<transbins_anu; n++)
      {
    chi2+=(cross_results_anu [k+1][l+1]/1e-41 - minerva_anu_results  [l][k])*
    (cross_results_anu [m+1][n+1]/1e-41 - minerva_anu_results  [n][m])*
    cov[l*longbins_anu+k  ] [n*longbins_anu+m ];
      }
    }
  }
}
nuwro<<"overall chi2 = ";
nuwro<<chi2<<endl;
cout<<"overall chi2 = "<<chi2<<endl;

//// now "small" chi2 in cosine bins 
  
/// a main loop (for chi2 only) in cosine


nuwro<<"small chi2 in cosine bins"<<endl;
  for (int jj=0; jj<longbins_anu; jj++)
  {
  
  for (int s=0; s<transbins_anu; s++)
  {
   for (int t=0; t<transbins_anu; t++)
   {
     cov1[s][t]=cov_copy [s*longbins_anu+jj] [t*longbins_anu+jj];//it is already invertible
   }
  }
cov1.Invert();
  
for (int k=0; k<transbins_anu; k++)
{
    for (int n=0; n<transbins_anu; n++)
    {
    chikwa[jj]+=(cross_results_anu [jj+1][k+1]/1e-41 - minerva_anu_results  [k][jj])*
    (cross_results_anu [jj+1][n+1]/1e-41 - minerva_anu_results  [n][jj])*
    cov1 [k][n];
  }
}

  }//end jj loop
  
  for (int kk=0; kk<longbins_anu; kk++)
  {cout<<chikwa[kk]<<"  ";
  nuwro<<kk<<" "<<chikwa[kk]<<"  ";
  }
  nuwro<<endl;
  
  
nuwro<<"chi2 without correlations"<<endl;

double chi22=0;
double errors [longbins_anu][transbins_anu];

for (int j=0; j<longbins_anu; j++)
  {
    for(int k=0; k<transbins_anu; k++)
    {
      if ( !(minerva_anu_results[k][j]==0) )
      {
	errors[j][k] = sqrt( minerva_anu_results_covariance [j+k*longbins_anu][j+k*longbins_anu] )*1e-41;
	
      chi22+= (minerva_anu_results[k][j]*1e-41 - cross_results_anu [j+1][k+1])*
      (minerva_anu_results[k][j]*1e-41 - cross_results_anu [j+1][k+1]) /errors[j][k]/errors[j][k];
      
      nuwro<<j<<" "<<k<<" "<<(minerva_anu_results[k][j]*1e-41 - cross_results_anu [j+1][k+1])*
      (minerva_anu_results[k][j]*1e-41 - cross_results_anu [j+1][k+1])
      /errors[j][k]/errors[j][k]<<endl;
      
      cout<<j<<" "<<k<<" "<<cross_results_anu [j+1][k+1]<<" "<<minerva_anu_results[k][j]*1e-41<<" "<<errors[j][k]<<endl;
      }
      else
      {nuwro<<j<<" "<<k<<" "<<endl;
      cout<<j<<" "<<k<<endl;
      }
    }
  }
nuwro<<endl;
nuwro<<"overall chi2 without correlations = "<<chi22<<endl;


///////////////////////////
/////// end chi2
//////////////////////////
  
double cross_mc=0;
double cross_data=0;
double sum_mc=0;
double sum_data=0;

	
///////////////////////////////////////////////
////////////////  Figures and normalizations/////////////////////
///////////////////////////////////////////////

nuwro<<"normalizations in bins"<<endl;
for (int j=0; j<longbins_anu; j++)
{
  double aa=long_anu_left[j]/1000;
  double bb=long_anu_left[j+1]/1000;
  int cc = int (long_anu_left[j]);
  
  sprintf (napiss, "sztukas_%d",j);
  sprintf (napiss2, "MINERvA_anu_0pi_pLpT_%d_mecscaling.png",j);
  sprintf (napiss3, "CC 0#pi #bar{#nu_{#mu}} %.1f<p_{L}<%.1f [GeV], #chi^{2}=%.2f",aa,bb, chikwa[j]);
  
  TH1D *hb = new TH1D (napiss, napiss, transbins_anu, trans_anu_left);
  
for (int k=0; k<transbins_anu; k++)
{
  xdpt_anu[k]=(trans_anu_left[k+1]+trans_anu_left[k])/2.0;
  exdpt_anu[k]=(trans_anu_left[k+1]-trans_anu_left[k])/2.0;
  ydpt_anu[k]=minerva_anu_results[k][j]*1e-41;
  eydpt_anu[k]=sqrt( minerva_anu_results_covariance [j+k*longbins_anu][j+k*longbins_anu] )*1e-41;
     hb->SetBinContent(k+1,cross_results_anu[j+1][k+1]);
}

for (int k=0; k<transbins_anu; k++)
{
  cross_mc   += cross_results_anu [j+1][k+1]   *(long_anu_left[j+1]-long_anu_left[j])*(trans_anu_left[k+1]-trans_anu_left[k])/1e6;
  cross_data += minerva_anu_results[k][j]*1e-41*(long_anu_left[j+1]-long_anu_left[j])*(trans_anu_left[k+1]-trans_anu_left[k])/1e6;
}

cout<<cross_mc<<"  "<<cross_data<<endl;
nuwro<<cross_mc<<"  "<<cross_data<<endl;
sum_mc+=cross_mc;
sum_data+=cross_data;
cross_mc=0;
cross_data=0;

TGraphErrors *gdpt_anu = new TGraphErrors(transbins_anu,xdpt_anu,ydpt_anu,exdpt_anu,eydpt_anu);

TCanvas *calphatb = new TCanvas("calphatb","calphatb",200,10,700,500);

   calphatb->Range(-0.2509977,-0.1717173,1.254176,1.076966);
   calphatb->SetFillColor(0);
   calphatb->SetBorderSize(2);
   calphatb->SetTickx();
   calphatb->SetTicky();
   calphatb->SetLeftMargin(0.1203828);//originally 0.13
   calphatb->SetRightMargin(0.02937799);//originally 0.06
   calphatb->SetTopMargin(0.0812855);//originally 0.06
   calphatb->SetBottomMargin(0.1375187);
   calphatb->SetFrameBorderMode(0);
   calphatb->SetFrameBorderMode(0);
   //calphatb->SetLogx();
   
hb->SetTitle(napiss3);
//hb->SetName("hb");
hb->SetLineColor(4);
hb->SetLineWidth(4);
hb->SetLineStyle(1);
hb->SetStats(kFALSE);
hb->SetName("nuwro");
hb->GetXaxis()->SetTitleSize(0.06);
hb->GetYaxis()->SetTitleSize(0.05);
hb->GetXaxis()->SetTitle("muon transverse momentum [MeV]");
hb->GetYaxis()->SetTitle("cross section (per nucleon) [cm^{2}/GeV^{2}/nucleon]");
//hb->GetYaxis()->SetRangeUser(0,3.8e-39);


if (j==0)
hb->GetYaxis()->SetRangeUser(0,2.1e-39);

if (j==1)
hb->GetYaxis()->SetRangeUser(0,2.6e-39);

if (j==2)
hb->GetYaxis()->SetRangeUser(0,2.75e-39);

if (j==3)
hb->GetYaxis()->SetRangeUser(0,2.4e-39);

if (j==4)
hb->GetYaxis()->SetRangeUser(0,1.4e-39);

if (j==5)
hb->GetYaxis()->SetRangeUser(0,0.7e-39);

if (j==6)
hb->GetYaxis()->SetRangeUser(0,0.3e-39);

if (j==7)
hb->GetYaxis()->SetRangeUser(0,0.15e-39);

if (j==8)
hb->GetYaxis()->SetRangeUser(0,75e-42);

if (j==9)
hb->GetYaxis()->SetRangeUser(0,37.5e-42);

hb->Draw();

gdpt_anu->SetName("minerva");
gdpt_anu->SetLineColor(kRed);
gdpt_anu->SetLineWidth(4);
gdpt_anu->Draw("E SAME"); 

TLegend *legendb = new TLegend(0.65,0.8,0.9,0.9);
   legendb->AddEntry("minerva","MINERvA", "lep");
   legendb->AddEntry("nuwro","NuWro", "L");
   legendb->Draw();

calphatb->Print(napiss2);
delete hb, calphatb, gdpt_anu, legendb;

}

nuwro<<"the final (overall) normalization check (first nuwro, then data)"<<endl;
cout<<sum_mc<<"  "<<sum_data<<endl;
nuwro<<sum_mc<<"  "<<sum_data<<endl;

//////////////////////////
///////// stack histograms
///////////////////////////

for (int j=0; j<longbins_anu; j++)
{
   double aa=long_anu_left[j]/1000;
  double bb=long_anu_left[j+1]/1000;
  int cc = int (long_anu_left[j]);

  sprintf (napiss, "sztukas_%d",j);
  sprintf (napiss2, "MINERvA_anu_0pi_pLpT_stack_%d_mecscaling.png",j);
  sprintf (napiss3, "CC 0#pi #bar{#nu_{#mu}}, stack, %.1f<p_{L}<%.1f [GeV], #chi^{2}=%.2f",aa,bb, chikwa[j]);
  
THStack *hs = new THStack("hs","Stacked 1D histograms");
hs->SetMinimum(0.0);

if (j==0)
hs->SetMaximum(2.1e-39);

if (j==1)
hs->SetMaximum(2.6e-39);

if (j==2)
hs->SetMaximum(2.75e-39);

if (j==3)
hs->SetMaximum(2.4e-39);

if (j==4)
hs->SetMaximum(1.4e-39);

if (j==5)
hs->SetMaximum(0.7e-39);

if (j==6)
hs->SetMaximum(0.3e-39);

if (j==7)
hs->SetMaximum(0.15e-39);

if (j==8)
hs->SetMaximum(75e-42);

if (j==9)
hs->SetMaximum(37.5e-42);

 TH1D *haqel = new TH1D ("haqel", "haqel", transbins_anu, trans_anu_left);
 TH1D *hares = new TH1D ("hares", "hares", transbins_anu, trans_anu_left);
 TH1D *hamec = new TH1D ("hamec", "hamec", transbins_anu, trans_anu_left);
 TH1D *hadis = new TH1D ("hadis", "hadis", transbins_anu, trans_anu_left);
 
   for (int k=0; k<longbins_anu; k++)
  {
    haqel->SetBinContent(k+1,cross_results_anu_qel[j+1][k+1]);
    hares->SetBinContent(k+1,cross_results_anu_res[j+1][k+1]);
    hamec->SetBinContent(k+1,cross_results_anu_mec[j+1][k+1]);
    hadis->SetBinContent(k+1,cross_results_anu_dis[j+1][k+1]);
  }

haqel->SetFillColor(kRed);
haqel->SetMarkerStyle(21);
haqel->SetMarkerColor(kRed);
haqel->SetName("ccqe");

hs->Add(haqel);

hares->SetFillColor(kBlue);
hares->SetMarkerStyle(21);
hares->SetMarkerColor(kBlue);
hares->SetName("res");

hs->Add(hares);

hamec->SetFillColor(kYellow);
hamec->SetMarkerStyle(21);
hamec->SetMarkerColor(kYellow);
hamec->SetName("mec");

hs->Add(hamec);

hadis->SetFillColor(kGreen);
hadis->SetMarkerStyle(21);
hadis->SetMarkerColor(kGreen);
hadis->SetName("dis");

hs->Add(hadis);

for (int k=0; k<transbins_anu; k++)
{
  xdpt_anu[k]=(trans_anu_left[k+1]+trans_anu_left[k])/2.0;
  exdpt_anu[k]=(trans_anu_left[k+1]-trans_anu_left[k])/2.0;
  ydpt_anu[k]=minerva_anu_results[k][j]*1e-41;
  eydpt_anu[k]=sqrt( minerva_anu_results_covariance [j+k*longbins_anu][j+k*longbins_anu] )*1e-41;
}

TGraphErrors *gpLpT = new TGraphErrors(transbins_anu,xdpt_anu,ydpt_anu,exdpt_anu,eydpt_anu);

TCanvas *calphat = new TCanvas("calphat","calphat",200,10,700,500);

   calphat->Range(-0.2509977,-0.1717173,1.254176,1.076966);
   calphat->SetFillColor(0);
   calphat->SetBorderSize(2);
   calphat->SetTickx();
   calphat->SetTicky();
   calphat->SetLeftMargin(0.1203828);//originally 0.13
   calphat->SetRightMargin(0.02937799);//originally 0.06
   calphat->SetTopMargin(0.0812855);//originally 0.06
   calphat->SetBottomMargin(0.1375187);
   calphat->SetFrameBorderMode(0);
   calphat->SetFrameBorderMode(0);

hs->Draw();
hs->SetTitle(napiss3);
hs->SetName("phit");
hs->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
hs->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
hs->GetXaxis()->SetTitle("muon transverse momentum [MeV]");
hs->GetYaxis()->SetTitle("cross section (per nucleon) [cm^{2}/GeV^{2}/nucleon]");
hs->Draw();

gpLpT->SetName("t2k");
gpLpT->SetLineWidth(3);
gpLpT->SetLineColor(kBlack);
gpLpT->Draw("E SAME"); 

TLegend *legendb = new TLegend(0.7,0.63,0.9,0.78);
   legendb->AddEntry("t2k","T2K", "lep");
   legendb->AddEntry("ccqe","NuWro CCQE", "f");
   legendb->AddEntry("res","NuWro RES", "f");
   legendb->AddEntry("mec","NuWro MEC", "f");
   legendb->AddEntry("dis","NuWro DIS", "f");
   legendb->Draw();
   
   calphat->Print(napiss2);
   delete hs, legendb, gpLpT, calphat;
}

/////////////////////////////////////////
///////// stack histograms complete mec
/////////////////////////////////////////

for (int j=0; j<longbins_anu; j++)
{
   double aa=long_anu_left[j]/1000;
  double bb=long_anu_left[j+1]/1000;
  int cc = int (long_anu_left[j]);

  sprintf (napiss, "sztukas_%d",j);
  sprintf (napiss2, "MINERvA_anu_0pi_pLpT_stack_%d_complete_mec_mecscaling.png",j);
  sprintf (napiss3, "CC 0#pi #bar{#nu_{#mu}}, stack, %.1f<p_{L}<%.1f [GeV], #chi^{2}=%.2f",aa,bb, chikwa[j]);
  
  sprintf (napiss4, "kinematical characteristics of MEC events (in bins from left energy and momentum transfer)");
  sprintf (napiss5, "1) w=%.0f#pm%.0f, 2) w=%.0f#pm%.0f, 3) w=%.0f#pm%.0f, 4) w=%.0f#pm%.0f, 5) w=%.0f#pm%.0f, 6) w=%.0f#pm%.0f",
	   w_average_mec[j][0],  w_sigma_mec[j][0], w_average_mec[j][1],  w_sigma_mec[j][1],
	   w_average_mec[j][2],  w_sigma_mec[j][2], w_average_mec[j][3],  w_sigma_mec[j][3], 
	   w_average_mec[j][4],  w_sigma_mec[j][4], w_average_mec[j][5],  w_sigma_mec[j][5]);
  sprintf (napiss6, "1) q=%.0f#pm%.0f, 2) q=%.0f#pm%.0f, 3) q=%.0f#pm%.0f, 4) q=%.0f#pm%.0f, 5) q=%.0f#pm%.0f, 6) q=%.0f#pm%.0f", 
	   q_average_mec[j][0],  q_sigma_mec[j][0], q_average_mec[j][1],  q_sigma_mec[j][1],
	   q_average_mec[j][2],  q_sigma_mec[j][2], q_average_mec[j][3],  q_sigma_mec[j][3], 
	   q_average_mec[j][4],  q_sigma_mec[j][4], q_average_mec[j][5],  q_sigma_mec[j][5]);
  
THStack *hs = new THStack("hs","Stacked 1D histograms");
hs->SetMinimum(0.0);

if (j==0)
hs->SetMaximum(2.45e-39);

if (j==1)
hs->SetMaximum(3e-39);

if (j==2)
hs->SetMaximum(3.15e-39);

if (j==3)
hs->SetMaximum(2.8e-39);

if (j==4)
hs->SetMaximum(1.6e-39);

if (j==5)
hs->SetMaximum(0.8e-39);

if (j==6)
hs->SetMaximum(0.35e-39);

if (j==7)
hs->SetMaximum(0.175e-39);

if (j==8)
hs->SetMaximum(86e-42);

if (j==9)
hs->SetMaximum(43e-42);

 TH1D *haqel = new TH1D ("haqel", "haqel", transbins_anu, trans_anu_left);
 TH1D *hares = new TH1D ("hares", "hares", transbins_anu, trans_anu_left);
 TH1D *hamec = new TH1D ("hamec", "hamec", transbins_anu, trans_anu_left);
 TH1D *hadis = new TH1D ("hadis", "hadis", transbins_anu, trans_anu_left);
 
   for (int k=0; k<longbins_anu; k++)
  {
    haqel->SetBinContent(k+1,cross_results_anu_qel[j+1][k+1]);
    hares->SetBinContent(k+1,cross_results_anu_res[j+1][k+1]);
    hamec->SetBinContent(k+1,cross_results_anu_mec[j+1][k+1]);
    hadis->SetBinContent(k+1,cross_results_anu_dis[j+1][k+1]);
  }

haqel->SetFillColor(kRed);
haqel->SetMarkerStyle(21);
haqel->SetMarkerColor(kRed);
haqel->SetName("ccqe");

hs->Add(haqel);

hares->SetFillColor(kBlue);
hares->SetMarkerStyle(21);
hares->SetMarkerColor(kBlue);
hares->SetName("res");

hs->Add(hares);

hamec->SetFillColor(kYellow);
hamec->SetMarkerStyle(21);
hamec->SetMarkerColor(kYellow);
hamec->SetName("mec");

hs->Add(hamec);

hadis->SetFillColor(kGreen);
hadis->SetMarkerStyle(21);
hadis->SetMarkerColor(kGreen);
hadis->SetName("dis");

hs->Add(hadis);

for (int k=0; k<transbins_anu; k++)
{
  xdpt_anu[k]=(trans_anu_left[k+1]+trans_anu_left[k])/2.0;
  exdpt_anu[k]=(trans_anu_left[k+1]-trans_anu_left[k])/2.0;
  ydpt_anu[k]=minerva_anu_results[k][j]*1e-41;
  eydpt_anu[k]=sqrt( minerva_anu_results_covariance [j+k*longbins_anu][j+k*longbins_anu] )*1e-41;
}

TGraphErrors *gpLpT = new TGraphErrors(transbins_anu,xdpt_anu,ydpt_anu,exdpt_anu,eydpt_anu);

TCanvas *calphat = new TCanvas("calphat","calphat",200,10,700,500);

   calphat->Range(-0.2509977,-0.1717173,1.254176,1.076966);
   calphat->SetFillColor(0);
   calphat->SetBorderSize(2);
   calphat->SetTickx();
   calphat->SetTicky();
   calphat->SetLeftMargin(0.1203828);//originally 0.13
   calphat->SetRightMargin(0.02937799);//originally 0.06
   calphat->SetTopMargin(0.0812855);//originally 0.06
   calphat->SetBottomMargin(0.1375187);
   calphat->SetFrameBorderMode(0);
   calphat->SetFrameBorderMode(0);

hs->Draw();
hs->SetTitle(napiss3);
hs->SetName("phit");
hs->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
hs->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
hs->GetXaxis()->SetTitle("muon transverse momentum [MeV]");
hs->GetYaxis()->SetTitle("cross section (per nucleon) [cm^{2}/GeV^{2}/nucleon]");
hs->Draw();

gpLpT->SetName("t2k");
gpLpT->SetLineWidth(3);
gpLpT->SetLineColor(kBlack);
gpLpT->Draw("E SAME"); 


TLegend *legendb = new TLegend(0.7,0.63,0.9,0.78);
   legendb->AddEntry("t2k","T2K", "lep");
   legendb->AddEntry("ccqe","NuWro CCQE", "f");
   legendb->AddEntry("res","NuWro RES", "f");
   legendb->AddEntry("mec","NuWro MEC", "f");
   legendb->AddEntry("dis","NuWro DIS", "f");
   legendb->Draw();
   
      TPaveText *pt = new TPaveText(0.15,0.79,0.95,0.91, "brNDC");
//   pt->AddText(napis5);
  pt->AddText(napiss4);
  pt->AddText(napiss5);
  pt->AddText(napiss6);
  //pt->AddText(napis8);
   pt->Draw();
   
   calphat->Print(napiss2);
   delete hs, legendb, gpLpT, calphat, pt;
}


/////////////////////////////////////////
///////// stack histograms complete qel
/////////////////////////////////////////

for (int j=0; j<longbins_anu; j++)
{
   double aa=long_anu_left[j]/1000;
  double bb=long_anu_left[j+1]/1000;
  int cc = int (long_anu_left[j]);

  sprintf (napiss, "sztukas_%d",j);
  sprintf (napiss2, "MINERvA_anu_0pi_pLpT_stack_%d_complete_qel_mecscaling.png",j);
  sprintf (napiss3, "CC 0#pi #bar{#nu_{#mu}}, stack, %.1f<p_{L}<%.1f [GeV], #chi^{2}=%.2f",aa,bb, chikwa[j]);
  
  sprintf (napisss4, "kinematical characteristics of CCQE events (in bins from left energy and momentum transfer)");
  sprintf (napisss5, "1) w=%.0f#pm%.0f, 2) w=%.0f#pm%.0f, 3) w=%.0f#pm%.0f, 4) w=%.0f#pm%.0f, 5) w=%.0f#pm%.0f, 6) w=%.0f#pm%.0f",
	   w_average_qel[j][0],  w_sigma_qel[j][0], w_average_qel[j][1],  w_sigma_qel[j][1],
	   w_average_qel[j][2],  w_sigma_qel[j][2], w_average_qel[j][3],  w_sigma_qel[j][3], 
	   w_average_qel[j][4],  w_sigma_qel[j][4], w_average_qel[j][5],  w_sigma_qel[j][5]);
  sprintf (napisss6, "1) q=%.0f#pm%.0f, 2) q=%.0f#pm%.0f, 3) q=%.0f#pm%.0f, 4) q=%.0f#pm%.0f, 5) q=%.0f#pm%.0f, 6) q=%.0f#pm%.0f", 
	   q_average_qel[j][0],  q_sigma_qel[j][0], q_average_qel[j][1],  q_sigma_qel[j][1],
	   q_average_qel[j][2],  q_sigma_qel[j][2], q_average_qel[j][3],  q_sigma_qel[j][3], 
	   q_average_qel[j][4],  q_sigma_qel[j][4], q_average_qel[j][5],  q_sigma_qel[j][5]);
  
THStack *hs = new THStack("hs","Stacked 1D histograms");
hs->SetMinimum(0.0);

if (j==0)
hs->SetMaximum(2.45e-39);

if (j==1)
hs->SetMaximum(3e-39);

if (j==2)
hs->SetMaximum(3.15e-39);

if (j==3)
hs->SetMaximum(2.8e-39);

if (j==4)
hs->SetMaximum(1.6e-39);

if (j==5)
hs->SetMaximum(0.8e-39);

if (j==6)
hs->SetMaximum(0.35e-39);

if (j==7)
hs->SetMaximum(0.175e-39);

if (j==8)
hs->SetMaximum(86e-42);

if (j==9)
hs->SetMaximum(43e-42);

 TH1D *haqel = new TH1D ("haqel", "haqel", transbins_anu, trans_anu_left);
 TH1D *hares = new TH1D ("hares", "hares", transbins_anu, trans_anu_left);
 TH1D *hamec = new TH1D ("hamec", "hamec", transbins_anu, trans_anu_left);
 TH1D *hadis = new TH1D ("hadis", "hadis", transbins_anu, trans_anu_left);
 
   for (int k=0; k<longbins_anu; k++)
  {
    haqel->SetBinContent(k+1,cross_results_anu_qel[j+1][k+1]);
    hares->SetBinContent(k+1,cross_results_anu_res[j+1][k+1]);
    hamec->SetBinContent(k+1,cross_results_anu_mec[j+1][k+1]);
    hadis->SetBinContent(k+1,cross_results_anu_dis[j+1][k+1]);
  }

haqel->SetFillColor(kRed);
haqel->SetMarkerStyle(21);
haqel->SetMarkerColor(kRed);
haqel->SetName("ccqe");

hs->Add(haqel);

hares->SetFillColor(kBlue);
hares->SetMarkerStyle(21);
hares->SetMarkerColor(kBlue);
hares->SetName("res");

hs->Add(hares);

hamec->SetFillColor(kYellow);
hamec->SetMarkerStyle(21);
hamec->SetMarkerColor(kYellow);
hamec->SetName("mec");

hs->Add(hamec);

hadis->SetFillColor(kGreen);
hadis->SetMarkerStyle(21);
hadis->SetMarkerColor(kGreen);
hadis->SetName("dis");

hs->Add(hadis);

for (int k=0; k<transbins_anu; k++)
{
  xdpt_anu[k]=(trans_anu_left[k+1]+trans_anu_left[k])/2.0;
  exdpt_anu[k]=(trans_anu_left[k+1]-trans_anu_left[k])/2.0;
  ydpt_anu[k]=minerva_anu_results[k][j]*1e-41;
  eydpt_anu[k]=sqrt( minerva_anu_results_covariance [j+k*longbins_anu][j+k*longbins_anu] )*1e-41;
}

TGraphErrors *gpLpT = new TGraphErrors(transbins_anu,xdpt_anu,ydpt_anu,exdpt_anu,eydpt_anu);

TCanvas *calphat = new TCanvas("calphat","calphat",200,10,700,500);

   calphat->Range(-0.2509977,-0.1717173,1.254176,1.076966);
   calphat->SetFillColor(0);
   calphat->SetBorderSize(2);
   calphat->SetTickx();
   calphat->SetTicky();
   calphat->SetLeftMargin(0.1203828);//originally 0.13
   calphat->SetRightMargin(0.02937799);//originally 0.06
   calphat->SetTopMargin(0.0812855);//originally 0.06
   calphat->SetBottomMargin(0.1375187);
   calphat->SetFrameBorderMode(0);
   calphat->SetFrameBorderMode(0);

hs->Draw();
hs->SetTitle(napiss3);
hs->SetName("phit");
hs->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
hs->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
hs->GetXaxis()->SetTitle("muon transverse momentum [MeV]");
hs->GetYaxis()->SetTitle("cross section (per nucleon) [cm^{2}/GeV^{2}/nucleon]");
hs->Draw();


gpLpT->SetName("t2k");
gpLpT->SetLineWidth(3);
gpLpT->SetLineColor(kBlack);
gpLpT->Draw("E SAME"); 


TLegend *legendb = new TLegend(0.7,0.63,0.9,0.78);
   legendb->AddEntry("t2k","T2K", "lep");
   legendb->AddEntry("ccqe","NuWro CCQE", "f");
   legendb->AddEntry("res","NuWro RES", "f");
   legendb->AddEntry("mec","NuWro MEC", "f");
   legendb->AddEntry("dis","NuWro DIS", "f");
   legendb->Draw();
   
TPaveText *ptpt = new TPaveText(0.15,0.79,0.95,0.91, "brNDC");
//   pt->AddText(napis5);
  ptpt->AddText(napisss4);
  ptpt->AddText(napisss5);
  ptpt->AddText(napisss6);
  //pt->AddText(napis8);
   ptpt->Draw();
   
   calphat->Print(napiss2);
   delete hs, legendb, gpLpT, calphat, ptpt;
}


}
