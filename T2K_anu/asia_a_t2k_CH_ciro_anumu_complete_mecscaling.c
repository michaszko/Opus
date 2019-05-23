#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

#include "TBrowser.h"
#include "TGraph2D.h"

//#include "t2k_CH_0pi_analysis1.h"
#include "ciro_data_new.h"
#include "mec_scaling_t2k_anu_april23_corr.h"
#include "mec_scaling_april23_corr.h"
#include "ciro_cov.h"

//to be changed
int events = 500000;

double Mproton = 938.272;
double Mproton2 = Mproton*Mproton;

double Mneutron = 939.565;

double mpiplus=139.57;
double mpiplus2=mpiplus*mpiplus;

double mmion=105.658;
double mmion2=mmion*mmion;

char napis [130];
char napis2 [130];
char napis3 [130];
char napis4 [130];
char napis5 [200];
char napis6 [200];
char napis7 [130];
char napis8 [130];

const int mombins = 14;
const int cosbins = 9;
double mom_left [mombins+1] = {5, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 2000, 3000, 5000, 30000};
double cos_left [cosbins+1] = {-1, 0.2, 0.6, 0.7, 0.8, 0.85, 0.9, 0.94, 0.98, 1};

double results [mombins+2] [cosbins+2];
double results_qel [mombins+2] [cosbins+2];
double results_res [mombins+2] [cosbins+2];
double results_mec [mombins+2] [cosbins+2];
double results_dis [mombins+2] [cosbins+2];

const int totalbins=58;
double cross [totalbins];
double cross_qel [totalbins];
double cross_res [totalbins];
double cross_mec [totalbins];
double cross_dis [totalbins];

int binning [cosbins] = {1, 5, 6, 6, 7, 8, 7, 10, 8};

int sklejka [totalbins] = {
  14,
  1,1,1,1,10,
  1,1,1,1,2,8,
  1,1,1,1,2,8,
  1,1,1,1,2,2,6,
  1,1,1,1,2,2,2,4,
  2,1,1,2,3,2,3,
  2,1,1,2,2,1,1,1,1,2,
  3,2,2,2,2,1,1,1
  };

  
int testsum1=0;
int testsum2=0;
  
double q_analyser_qel [mombins+2][cosbins+2];
double q2_analyser_qel [mombins+2][cosbins+2];
double q_analyser_mec [mombins+2][cosbins+2];
double q2_analyser_mec [mombins+2][cosbins+2];
double w_analyser_qel [mombins+2][cosbins+2];
double w2_analyser_qel [mombins+2][cosbins+2];
double w_analyser_mec [mombins+2][cosbins+2];
double w2_analyser_mec [mombins+2][cosbins+2];

double q_average_qel [totalbins];
double q_sigma_qel [totalbins];
double w_average_qel [totalbins];
double w_sigma_qel [totalbins];

double q_average_mec [totalbins];
double q_sigma_mec [totalbins];
double w_average_mec [totalbins];
double w_sigma_mec [totalbins];

///////////////////////////////////////////////////////
double mscaling (double momtran, double entran)
{
  int momnumber = int(momtran/100);
  int ennumber = int(entran/100);
  
  //cout<<momtran<<"  "<<momnumber<<"  "<<entran<<"  "<<ennumber<<"  "<<mec_scaling [momnumber][ennumber]<<endl;
  //return mec_scaling [momnumber][ennumber];
  //return mec_scaling_t2k_anu [momnumber][ennumber];
  return 1;
}
//////////////////////////////////////////////////////

ofstream nuwro ("1709_t2k_CH_ciro_anumu_report_oct5.dat");
//ofstream nuwro2 ("inverse_cov_anu.txt");

void asia_a (TFile *Input, int ilezdarzen, double table[mombins+2] [cosbins+2], 
	     double table_qel[mombins+2] [cosbins+2], double table_res[mombins+2] [cosbins+2], 
	     double table_mec[mombins+2] [cosbins+2], double table_dis[mombins+2] [cosbins+2], 
	     double table_qel_w [mombins+2] [cosbins+2], 
	     double table_qel_w2 [mombins+2] [cosbins+2], double table_qel_q [mombins+2] [cosbins+2], 
	     double table_qel_q2 [mombins+2] [cosbins+2], double table_mec_w [mombins+2] [cosbins+2], 
	     double table_mec_w2 [mombins+2] [cosbins+2], double table_mec_q [mombins+2] [cosbins+2], 
	     double table_mec_q2 [mombins+2] [cosbins+2])
{
        TH2D *h1 = new TH2D ("h1", "h1", mombins, mom_left, cosbins, cos_left); 
	TH2D *h1qel = new TH2D ("h1qel", "h1qel", mombins, mom_left, cosbins, cos_left); 
	TH2D *h1res = new TH2D ("h1res", "h1res", mombins, mom_left, cosbins, cos_left); 
	TH2D *h1mec = new TH2D ("h1mec", "h1mec", mombins, mom_left, cosbins, cos_left); 
	TH2D *h1dis = new TH2D ("h1dis", "h1dis", mombins, mom_left, cosbins, cos_left);
	
	TH2D *h1qel_w = new TH2D ("h1qel_w", "h1qel_w", mombins, mom_left, cosbins, cos_left); 
	TH2D *h1qel_w2 = new TH2D ("h1qel_w2", "h1qel_w2", mombins, mom_left, cosbins, cos_left); 
	TH2D *h1qel_q = new TH2D ("h1qel_q", "h1qel_q", mombins, mom_left, cosbins, cos_left); 
	TH2D *h1qel_q2 = new TH2D ("h1qel_q2", "h1qel_q2", mombins, mom_left, cosbins, cos_left); 
	
	TH2D *h1mec_w = new TH2D ("h1mec_w", "h1mec_w", mombins, mom_left, cosbins, cos_left); 
	TH2D *h1mec_w2 = new TH2D ("h1mec_w2", "h1mec_w2", mombins, mom_left, cosbins, cos_left); 
	TH2D *h1mec_q = new TH2D ("h1mec_q", "h1mec_q", mombins, mom_left, cosbins, cos_left); 
	TH2D *h1mec_q2 = new TH2D ("h1mec_q2", "h1mec_q2", mombins, mom_left, cosbins, cos_left); 
	
	TTree *tt1 = (TTree*)Input->Get("treeout");
	event * e= new event();
	tt1->SetBranchAddress("e",&e);
		
	for( int k=0; k < ilezdarzen; k++ )
	{
	tt1->GetEntry(k);
		{ 
		  if (e->fof(211, -211, 111)==0)//no pions 
		  {
		  double muoncos = e->out[0].z/e->out[0].momentum();
		  double muonmom = e->out[0].momentum();
		   double entrans = e->in[0].t - e->out[0].t;
		   double momtrans2 = entrans*entrans - e->q2();
		   double momtrans = sqrt(momtrans2);
		   
		   if (e->flag.mec==1)
		  h1->Fill(muonmom, muoncos,mscaling(momtrans,entrans) );
		  else
		    h1->Fill(muonmom, muoncos);
		  
		   if (e->flag.qel==1)
		      {h1qel->Fill(muonmom, muoncos);
			h1qel_w->Fill(muonmom, muoncos,entrans);
			h1qel_w2->Fill(muonmom, muoncos,entrans*entrans);
			h1qel_q->Fill(muonmom, muoncos,sqrt(momtrans2));
			h1qel_q2->Fill(muonmom, muoncos,momtrans2);
		      }
		       if (e->flag.res==1)
			h1res->Fill(muonmom, muoncos);
		        if (e->flag.mec==1)
			{h1mec->Fill(muonmom, muoncos,mscaling(momtrans,entrans));
			h1mec_w->Fill(muonmom, muoncos,entrans);
			h1mec_w2->Fill(muonmom, muoncos,entrans*entrans);
			h1mec_q->Fill(muonmom, muoncos,sqrt(momtrans2));
			h1mec_q2->Fill(muonmom, muoncos,momtrans2);
			}
		        if (e->flag.dis==1)
			h1dis->Fill(muonmom, muoncos);
		  /*
		  if (e->flag.qel==1)
		    h1qel->Fill(muonmom, muoncos);
		  if (e->flag.res==1)
		    h1res->Fill(muonmom, muoncos);
		  if (e->flag.mec==1)
		    h1mec->Fill(muonmom, muoncos);
		  if (e->flag.dis==1)
		    h1dis->Fill(muonmom, muoncos);
		  */
		  }
		}
	}
	delete e;
	
	for (int j=0; j<mombins+2; j++)
	{
	  for (int k=0; k<cosbins+2; k++)
	  { /*
	    table[j][k] = h1->GetBinContent(j,k);
	    table_qel[j][k] = h1qel->GetBinContent(j,k);
	    table_res[j][k] = h1res->GetBinContent(j,k);
	    table_mec[j][k] = h1mec->GetBinContent(j,k);
	    table_dis[j][k] = h1dis->GetBinContent(j,k);
	    */
	    table [j][k] = h1->GetBinContent(j,k);
	    table_qel [j][k] = h1qel->GetBinContent(j,k);
	    table_qel_w [j][k] = h1qel_w->GetBinContent(j,k);
	    table_qel_w2 [j][k] = h1qel_w2->GetBinContent(j,k);
	    table_qel_q [j][k] = h1qel_q->GetBinContent(j,k);
	    table_qel_q2 [j][k] = h1qel_q2->GetBinContent(j,k);
	    table_res [j][k] = h1res->GetBinContent(j,k);
	    table_mec [j][k] = h1mec->GetBinContent(j,k);
	    table_mec_w [j][k] = h1mec_w->GetBinContent(j,k);
	    table_mec_w2 [j][k] = h1mec_w2->GetBinContent(j,k);
	    table_mec_q [j][k] = h1mec_q->GetBinContent(j,k);
	    table_mec_q2 [j][k] = h1mec_q2->GetBinContent(j,k);
	    table_dis [j][k] = h1dis->GetBinContent(j,k);
	  }
	}
	delete h1, h1qel, h1res, h1mec, h1dis, h1qel_w, h1qel_w2, h1qel_q, h1qel_q2, h1mec_w, h1mec_w2, h1mec_q, h1mec_q2;
}

double norm_cc (ifstream &Input)
{
string ff;
double answer=0;int testsum1=0;
int tests
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
void asia_a_t2k_CH_ciro_anumu_complete_mecscaling()
/////////////////////////////////////////
////////////////////////////////////////
{
ifstream Input ("a_t2k_cc_CH_lfg_500kilo_anu.root.txt");
double cross_t2k = norm_cc(Input);

TFile *tf1 = new TFile("a_t2k_cc_CH_lfg_500kilo_anu.root");
asia_a (tf1, events, results, results_qel, results_res, results_mec, results_dis, w_analyser_qel, w2_analyser_qel, 
	q_analyser_qel, q2_analyser_qel, w_analyser_mec, w2_analyser_mec, 
	q_analyser_mec, q2_analyser_mec);
delete tf1;

int aggr [totalbins];//how many bins before

aggr[0]=0;
for (int j=1; j<totalbins; j++)
aggr[j]=aggr[j-1]+sklejka[j-1];

////////////////////////////////
//////// cross sections, errors
////////////////////////////////
  
nuwro<<"cross section results"<<endl;
for (int j=0; j<totalbins; j++)
{
    int dlugosc = sklejka [j];
    int koss = aggr[j]/mombins;//which cos bin
    int momm = aggr[j] - koss*mombins;//which mom bin
    double bin_start = mom_left[momm];
    double bin_end =mom_left[momm+dlugosc];//correct understanding of binning
  
    double sumka= 0;
    double sumka_qel= 0;
    double sumka_res= 0;
    double sumka_mec= 0;
    double sumka_dis= 0;
    
    double sumka_qel_q=0;
    double sumka_qel_q2=0;
    double sumka_qel_w=0;
    double sumka_qel_w2=0;
    
    double sumka_mec_q=0;
    double sumka_mec_q2=0;
    double sumka_mec_w=0;
    double sumka_mec_w2=0;
    
  for (int s=momm+1; s<momm+dlugosc+1; s++)
  {sumka += results[s][koss+1];
    sumka_qel += results_qel[s][koss+1];
    sumka_res += results_res[s][koss+1];
    sumka_mec += results_mec[s][koss+1];
    sumka_dis += results_dis[s][koss+1];
    
    sumka_qel_q += q_analyser_qel[s][koss+1];
    sumka_qel_q2 += q2_analyser_qel[s][koss+1];
    sumka_qel_w += w_analyser_qel[s][koss+1];
    sumka_qel_w2 += w2_analyser_qel[s][koss+1];
    sumka_mec_q += q_analyser_mec[s][koss+1];
    sumka_mec_q2 += q2_analyser_mec[s][koss+1];
    sumka_mec_w += w_analyser_mec[s][koss+1];
    sumka_mec_w2 += w2_analyser_mec[s][koss+1];
  }
  
  cross[j]=sumka/events*cross_t2k/(cos_left[koss+1]-cos_left[koss])/(mom_left[momm+dlugosc]-mom_left[momm])*1e3;
  cross_qel[j]=sumka_qel/events*cross_t2k/(cos_left[koss+1]-cos_left[koss])/(mom_left[momm+dlugosc]-mom_left[momm])*1e3;
  cross_res[j]=sumka_res/events*cross_t2k/(cos_left[koss+1]-cos_left[koss])/(mom_left[momm+dlugosc]-mom_left[momm])*1e3;
  cross_mec[j]=sumka_mec/events*cross_t2k/(cos_left[koss+1]-cos_left[koss])/(mom_left[momm+dlugosc]-mom_left[momm])*1e3;
  cross_dis[j]=sumka_dis/events*cross_t2k/(cos_left[koss+1]-cos_left[koss])/(mom_left[momm+dlugosc]-mom_left[momm])*1e3;
  //nuwro<<cross[j]<<"  "<<cross_qel[j]<<"  "<<cross_res[j]<<"  "<<cross_mec[j]<<"  "<<cross_dis[j]<<"  ";
  
    double srednia_w_qel;
    double dyspersja_w_qel;
    double srednia_q_qel;
    double dyspersja_q_qel;
    
    if (sumka_qel>0)
    {
      srednia_w_qel= sumka_qel_w/sumka_qel;
      dyspersja_w_qel = sqrt( sumka_qel_w2/sumka_qel - srednia_w_qel*srednia_w_qel );
      srednia_q_qel= sumka_qel_q/sumka_qel;
      dyspersja_q_qel = sqrt( sumka_qel_q2/sumka_qel - srednia_q_qel*srednia_q_qel );
    }
    else
    {
      srednia_w_qel =0;
      dyspersja_w_qel=0;
      srednia_q_qel =0;
      dyspersja_q_qel=0;
    }
q_average_qel[j]= srednia_q_qel;
q_sigma_qel[j]= dyspersja_q_qel;
w_average_qel[j]= srednia_w_qel;
w_sigma_qel[j]= dyspersja_w_qel;

    double srednia_w_mec;
    double dyspersja_w_mec;
    double srednia_q_mec;
    double dyspersja_q_mec; 
    
    if (sumka_mec>0)
    {
      srednia_w_mec= sumka_mec_w/sumka_mec;
      dyspersja_w_mec = sqrt( sumka_mec_w2/sumka_mec - srednia_w_mec*srednia_w_mec );
      srednia_q_mec= sumka_mec_q/sumka_mec;
      dyspersja_q_mec = sqrt( sumka_mec_q2/sumka_mec - srednia_q_mec*srednia_q_mec );
    }
    else
    {
      srednia_w_mec =0;
      dyspersja_w_mec=0;
      srednia_q_mec =0;
      dyspersja_q_mec=0;
    }
q_average_mec[j]= srednia_q_mec;
q_sigma_mec[j]= dyspersja_q_mec;
w_average_mec[j]= srednia_w_mec;
w_sigma_mec[j]= dyspersja_w_mec;
    
  }
  
nuwro<<"kinematical analysis"<<endl;

for (int t=0; t<totalbins; t++)
  nuwro<<q_average_qel[t]<<"  "<<q_sigma_qel[t]<<"  "<<w_average_qel[t]<<"  "<<w_sigma_qel[t]<<endl;
nuwro<<"end of kinematical analysis"<<endl;

double error[totalbins];
nuwro<<"cross sections and errors"<<endl;

for (int s=0; s<totalbins; s++)
{error [s] = //sqrt( //ciro_anu_cov_norm [s*(totalbins+1)]+ ciro_anu_cov_shape [s*(totalbins+1)] + 
  ciro_anu_errors[s];
nuwro<<cross[s]<<"  "<<ciro_anu_sections[s]<<"  "<<error[s]<<endl;
}
nuwro<<endl;

///////////////////////////////////////////////
////////////////  chi2 /////////////////////
///////////////////////////////////////////////

// overall chi2 without correlations

double chi2=0;


for (int j=0; j<totalbins; j++)
{
  chi2 += (ciro_anu_sections [j]-cross[j])*(ciro_anu_sections [j]-cross[j])/ciro_anu_errors[j]/ciro_anu_errors[j];
  nuwro<<j<<"  "<<(ciro_anu_sections [j]-cross[j])*(ciro_anu_sections [j]-cross[j])/ciro_anu_errors[j]/ciro_anu_errors[j]<<endl;
}

nuwro<<endl;
nuwro<<"overall chi2 neglecting correlations = "<<chi2<<endl;
nuwro<<endl;


 
int sumcosines=0;
int cosinewidth=0;

TMatrixD cov(totalbins,totalbins);
TMatrixD cov_copy(totalbins,totalbins);
double chi2cov=0;
double wklad;

for (int k=0; k<totalbins; k++)
  {
    for (int j=0; j<totalbins; j++)
    {
      cov [k][j] = //ciro_anu_cov_norm [k*totalbins+j] + ciro_anu_cov_shape [k*totalbins+j] + 
      ciro_cov [58+k][58+j];
      cov_copy[k][j] = cov[k][j];
    }
  }

cov.Invert();

/*
for (int k=0; k<totalbins; k++)
  {
    for (int j=0; j<totalbins; j++)
    {
      nuwro2<<cov[k][j]<<",  ";
    }
    nuwro2<<endl;
  }
  nuwro2<<endl;
*/


nuwro<<"chi2 study"<<endl;
for (int k=0; k<totalbins; k++)
{
  for(int l=0; l<totalbins; l++)
  {
    wklad=( cross [k]/1e-39 - ciro_anu_sections [k]/1e-39 )*( cross [l]/1e-39- ciro_anu_sections[l]/1e-39 )*
    cov[k][l];
    chi2cov+=wklad;
    //if (wklad >100)
      //nuwro<<k<<"  "<<l<<"  "<<wklad<<endl;
  }
}
nuwro<<"overall chi2 computation = "<<chi2cov<<endl;

/*
int used=0;
nuwro<<"small chi2 computations"<<endl;
double chikwa[cosbins];

for (int n=0; n<cosbins; n++)
  chikwa[n]=0;//starting point for chi2 in cosine bins
 
for (int s=0; s<cosbins; s++)//computation of "small" chi2
{
  int wymiar = binning[s];
  TMatrixD cov2(wymiar,wymiar);
  
  for (int m=0; m<wymiar; m++)
  {
    for (int n=0; n<wymiar; n++)
    {
      cov2[m][n]=cov_copy[used+m][used+n];
    }
  }
  
cov2.Invert();

for (int m=0; m<wymiar; m++)
  {
    for (int n=0; n<wymiar; n++)
    {
      chikwa[s]+=(cross [used+m]/1e-39 - ciro_anu_sections [used+m]/1e-39)*
      (cross [used+n]/1e-39 - ciro_anu_sections [used+n]/1e-39)*cov2[m][n];
    }
  }
nuwro<<chikwa[s]<<"  ";
used+=binning[s];
}
nuwro<<"end of chi2"<<endl;
*/
///////////////////////////////////////////////
////////////////  Figures /////////////////////
///////////////////////////////////////////////

int pointer=0;

double cross_mc=0;
double cross_data=0;
double sum_mc=0;
double sum_data=0;
double partialsum=0;

nuwro<<"normalizations (nuwro/data); the last bin is subtracted"<<endl;
for (int t=0; t<cosbins; t++)//figures in cosine bins
{
  
for (int j=0; j<totalbins; j++)//identification of the first momentum bin
{
 if (aggr[j]==mombins*t)
   pointer=j;
}
  
  double a=cos_left[t];
  double b=cos_left[t+1];//start end of the cos bin, this part is easy
  int c = int (100*a+100);
  
  sprintf (napis, "sztuka_%d",t);
  sprintf (napis2, "T2K_CC0pi_ciro_anumu_%d_oct5.png",c);
  //sprintf (napis3, "T2K CC0#pi CH, new, anumu, %.2f<cos#theta_{#mu}<%.2f, #chi^{2}=%.1f",a,b,chikwa[t]);
  sprintf (napis3, "T2K CC0#pi CH, new, anumu, %.2f<cos#theta_{#mu}<%.2f",a,b);
  //int mombins_new=binning[t]-1;//number of momentum bins
  int mombins_new=binning[t]-1;//number of momentum bins
  double mom_left_new [mombins_new+1];//array with new binning
  
  mom_left_new[0]=5;//to avoid problems in log scale plots
  int mom_run=0;
  
  for (int s=1; s<mombins_new+1; s++)
  {mom_left_new[s]=mom_left[sklejka[pointer]+mom_run];
    mom_run+=sklejka[pointer+s];
  }
  
double xdpt[mombins_new];
double ydpt[mombins_new];
double exdpt[mombins_new];
double eydpt[mombins_new];

for (int k=0; k<mombins_new; k++)
{
  cross_mc   += cross[pointer+k]     *(mom_left_new[k+1]-mom_left_new[k])*(cos_left[t+1]-cos_left[t])/1e3;
  cross_data += ciro_anu_sections[pointer+k]*(mom_left_new[k+1]-mom_left_new[k])*(cos_left[t+1]-cos_left[t])/1e3;
}

nuwro<<t<<"  "<<cross_mc<<"  "<<cross_data<<endl;
sum_mc+=cross_mc;
sum_data+=cross_data;
cross_mc=0;
cross_data=0;

  TH1D *ha = new TH1D (napis, napis, mombins_new, mom_left_new);
  
   for (int k=0; k<mombins_new; k++)
   {ha->SetBinContent(k+1,cross[pointer+k]);
     
    xdpt[k]=(mom_left_new[k+1]+mom_left_new[k])/2.0;
    exdpt[k]=(mom_left_new[k+1]-mom_left_new[k])/2.0;
    ydpt[k]= ciro_anu_sections[pointer+k];
    eydpt[k]=error[pointer+k];
   }
   
   TGraphErrors *gdpt = new TGraphErrors(mombins_new,xdpt,ydpt,exdpt,eydpt);
    
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
   calphat->SetLogx();

ha->SetTitle(napis3);
ha->SetName("ha");
ha->SetLineColor(4);
ha->SetLineWidth(4);
ha->SetLineStyle(1);
ha->SetStats(kFALSE);
ha->SetName("nuwro");
ha->GetXaxis()->SetTitleSize(0.06);
ha->GetYaxis()->SetTitleSize(0.05);
ha->GetXaxis()->SetTitle("muon momentum [MeV]");
ha->GetYaxis()->SetTitle("cross section [cm^{2}/GeV]");

if (t==0)
ha->GetYaxis()->SetRangeUser(0,3.2e-39);

if (t==1)
ha->GetYaxis()->SetRangeUser(0,1.6e-39);

if (t==2)
ha->GetYaxis()->SetRangeUser(0,5e-39);

if (t==3)
ha->GetYaxis()->SetRangeUser(0,5e-39);

if (t==4)
ha->GetYaxis()->SetRangeUser(0,6.5e-39);

if (t==5)
ha->GetYaxis()->SetRangeUser(0,7e-39);

if (t==6)
ha->GetYaxis()->SetRangeUser(0,8e-39);
 
if (t==7)
ha->GetYaxis()->SetRangeUser(0,8e-39);
  
if (t==8)
ha->GetYaxis()->SetRangeUser(0,5e-39);
  
ha->Draw();

gdpt->SetName("t2k");
gdpt->SetLineColor(kViolet);
gdpt->SetLineWidth(4);
gdpt->Draw("E SAME"); 

/*if  (t>3)
{TLegend *legenda = new TLegend(0.65,0.8,0.9,0.9);
legenda->AddEntry("t2k","T2K", "lep");
   legenda->AddEntry("nuwro","NuWro", "L");
   legenda->Draw();
}*/
if (t<9)
{TLegend *legenda = new TLegend(0.15,0.8,0.4,0.9);
  
   legenda->AddEntry("t2k","T2K", "lep");
   legenda->AddEntry("nuwro","NuWro", "L");
   legenda->Draw();
}
calphat->Print(napis2);
delete ha, calphat, gdpt;

}
nuwro<<endl;
nuwro<<"overall normalization"<<endl;
nuwro<<sum_mc<<"  "<<sum_data<<endl;

/////////////////////////////////////
///////////// stack //////////////////
/////////////////////////////////////

pointer=0;

for (int t=0; t<cosbins; t++)//figures in cosine bins
{
  
for (int j=0; j<totalbins; j++)//identification of the first momentum bin
{
 if (aggr[j]==mombins*t)
   pointer=j;
}

  double a=cos_left[t];
  double b=cos_left[t+1];//start end of the cos bin, this part is easy
  int c = int (100*a+100);
  
  sprintf (napis, "sztuka_%d",t);
  sprintf (napis2, "T2K_CC0pi_ciro_anumu_%d_stack_oct5.png",c);
  //sprintf (napis3, "T2K CC0#pi CH, new, anumu, %.2f<cos#theta_{#mu}<%.2f, #chi^{2}=%.1f",a,b,chikwa[t]);
  sprintf (napis3, "T2K CC0#pi CH, new, anumu, %.2f<cos#theta_{#mu}<%.2f",a,b);
  
  int mombins_new=binning[t]-1;//number of momentum bins
  double mom_left_new [mombins_new+1];//array with new binning
  
  mom_left_new[0]=5;//to avoid problems in log scale plots
  int mom_run=0;
  
  for (int s=1; s<mombins_new+1; s++)
  {mom_left_new[s]=mom_left[sklejka[pointer]+mom_run];
    mom_run+=sklejka[pointer+s];
  }
 
THStack *hs = new THStack("hs","Stacked 1D histograms");
hs->SetMinimum(0.0);

if (t==0)
hs->SetMaximum(3.2e-39);

if (t==1)
hs->SetMaximum(1.6e-39);

if (t==2)
hs->SetMaximum(5e-39);

if (t==3)
hs->SetMaximum(5e-39);

if (t==4)
hs->SetMaximum(6.5e-39);

if (t==5)
hs->SetMaximum(7e-39);

if (t==6)
hs->SetMaximum(8e-39);
 
if (t==7)
hs->SetMaximum(8e-39);
  
if (t==8)
hs->SetMaximum(5e-39);
   
TH1D *haqel = new TH1D ("haqel", "haqel", mombins_new, mom_left_new);
 TH1D *hares = new TH1D ("hares", "hares", mombins_new, mom_left_new);
 TH1D *hamec = new TH1D ("hamec", "hamec", mombins_new, mom_left_new);
 TH1D *hadis = new TH1D ("hadis", "hadis", mombins_new, mom_left_new);
 
double xdpt[mombins_new];
double ydpt[mombins_new];
double exdpt[mombins_new];
double eydpt[mombins_new];
 
for (int k=0; k<mombins_new; k++)
   {haqel->SetBinContent(k+1,cross_qel[pointer+k]);
     hares->SetBinContent(k+1,cross_res[pointer+k]);
     hamec->SetBinContent(k+1,cross_mec[pointer+k]);
     hadis->SetBinContent(k+1,cross_dis[pointer+k]);
     
    xdpt[k]=(mom_left_new[k+1]+mom_left_new[k])/2.0;
    exdpt[k]=(mom_left_new[k+1]-mom_left_new[k])/2.0;
    ydpt[k]= ciro_anu_sections[pointer+k];
    eydpt[k]=error[pointer+k];
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

TGraphErrors *gan2 = new TGraphErrors(mombins_new,xdpt,ydpt,exdpt,eydpt);
 

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
   calphat->SetLogx();

hs->Draw();
hs->SetTitle(napis3);
hs->SetName("phit");
hs->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
hs->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
hs->GetXaxis()->SetTitle("muon momentum [MeV]");
hs->GetYaxis()->SetTitle("cross section (per nucleon) [cm^{2}/GeV^{2}/nucleon]");
hs->Draw();

gan2->SetName("t2k");
gan2->SetLineWidth(3);
gan2->SetLineColor(kBlack);
gan2->Draw("E SAME"); 

if (t<9)
{TLegend *legendb = new TLegend(0.15,0.7,0.4,0.9);
   legendb->AddEntry("t2k","T2K", "lep");
   legendb->AddEntry("ccqe","NuWro CCQE", "f");
   legendb->AddEntry("res","NuWro RES", "f");
   legendb->AddEntry("mec","NuWro MEC", "f");
   legendb->AddEntry("dis","NuWro DIS", "f");
   legendb->Draw();
}
/*
if (t>4)
{TLegend *legendb = new TLegend(0.65,0.7,0.9,0.9);
   legendb->AddEntry("t2k","T2K", "lep");
   legendb->AddEntry("ccqe","NuWro CCQE", "f");
   legendb->AddEntry("res","NuWro RES", "f");
   legendb->AddEntry("mec","NuWro MEC", "f");
   legendb->AddEntry("dis","NuWro DIS", "f");
   legendb->Draw();
}
  */
   calphat->Print(napis2);
   delete hs, gan2, calphat, haqel, hares, hamec, hadis;
}


//////////////////////////
///////// stack (complete) histograms (qel)
///////////////////////////

pointer=0;
int used2=0;
for (int t=0; t<cosbins; t++)//figures in cosine bins
{
  
for (int j=0; j<totalbins; j++)//identification of the first momentum bin
{
 if (aggr[j]==mombins*t)
   pointer=j;
}

  double a=cos_left[t];
  double b=cos_left[t+1];//start end of the cos bin, this part is easy
  int c = int (100*a+100);
  
  int mombins_new=binning[t]-1;//number of momentum bins
  double mom_left_new [mombins_new+1];//array with new binning
  
  mom_left_new[0]=5;//to avoid problems in log scale plots
  int mom_run=0;
  
  for (int s=1; s<mombins_new+1; s++)
  {mom_left_new[s]=mom_left[sklejka[pointer]+mom_run];
    mom_run+=sklejka[pointer+s];
  }
  
sprintf (napis, "sztuka_%d",t);
  sprintf (napis2, "T2K_CC0pi_ciro_anumu_%d_stack_complete_qel.png",c);
  //sprintf (napis3, "T2K CC0#pi CH, new, anumu, %.2f<cos#theta_{#mu}<%.2f, #chi^{2}=%.1f",a,b,chikwa[t]);
  sprintf (napis3, "T2K CC0#pi CH, new, anumu, %.2f<cos#theta_{#mu}<%.2f",a,b);
  sprintf (napis4, "kinematical characteristics of CCQE events (in bins from left energy and momentum transfer)");
  
  sprintf (napis5, "1) w=%.0f#pm%.0f, ",w_average_qel[used2], w_sigma_qel[used2]);
  for (int s=1; s<mombins_new; s++)
sprintf(napis5 + strlen(napis5),"%d) w=%.0f#pm%.0f, ",s+1, w_average_qel[used2+s], w_sigma_qel[used2+s]);

 sprintf (napis6, "1) q=%.0f#pm%.0f, ",q_average_qel[used2], q_sigma_qel[used2]);
  for (int t=1; t<mombins_new; t++)
sprintf(napis6 + strlen(napis6),"%d) q=%.0f#pm%.0f, ",t+1, q_average_qel[used2+t], q_sigma_qel[used2+t]);

THStack *hs = new THStack("hs","Stacked 1D histograms");
hs->SetMinimum(0.0);

if (t==0)
hs->SetMaximum(4.9e-39);

if (t==1)
hs->SetMaximum(2e-39);

if (t==2)
hs->SetMaximum(6.5e-39);

if (t==3)
hs->SetMaximum(6.5e-39);

if (t==4)
hs->SetMaximum(8e-39);

if (t==5)
hs->SetMaximum(9e-39);

if (t==6)
hs->SetMaximum(10e-39);
 
if (t==7)
hs->SetMaximum(10e-39);
  
if (t==8)
hs->SetMaximum(6.25e-39);
   
TH1D *haqel = new TH1D ("haqel", "haqel", mombins_new, mom_left_new);
 TH1D *hares = new TH1D ("hares", "hares", mombins_new, mom_left_new);
 TH1D *hamec = new TH1D ("hamec", "hamec", mombins_new, mom_left_new);
 TH1D *hadis = new TH1D ("hadis", "hadis", mombins_new, mom_left_new);
 
double xdpt[mombins_new];
double ydpt[mombins_new];
double exdpt[mombins_new];
double eydpt[mombins_new];
 
for (int k=0; k<mombins_new; k++)
   {haqel->SetBinContent(k+1,cross_qel[pointer+k]);
     hares->SetBinContent(k+1,cross_res[pointer+k]);
     hamec->SetBinContent(k+1,cross_mec[pointer+k]);
     hadis->SetBinContent(k+1,cross_dis[pointer+k]);
     
    xdpt[k]=(mom_left_new[k+1]+mom_left_new[k])/2.0;
    exdpt[k]=(mom_left_new[k+1]-mom_left_new[k])/2.0;
    ydpt[k]= ciro_anu_sections[pointer+k];
    eydpt[k]=error[pointer+k];
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

TGraphErrors *gan2 = new TGraphErrors(mombins_new,xdpt,ydpt,exdpt,eydpt);
 

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
   calphat->SetLogx();

hs->Draw();
hs->SetTitle(napis3);
hs->SetName("phit");
hs->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
hs->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
hs->GetXaxis()->SetTitle("muon momentum [MeV]");
hs->GetYaxis()->SetTitle("cross section (per nucleon) [cm^{2}/GeV^{2}/nucleon]");
hs->Draw();

gan2->SetName("t2k");
gan2->SetLineWidth(3);
gan2->SetLineColor(kBlack);
gan2->Draw("E SAME"); 

if (t<5)
{TLegend *legendb = new TLegend(0.15,0.6,0.4,0.75);
   legendb->AddEntry("t2k","T2K", "lep");
   legendb->AddEntry("ccqe","NuWro CCQE", "f");
   legendb->AddEntry("res","NuWro RES", "f");
   legendb->AddEntry("mec","NuWro MEC", "f");
   legendb->AddEntry("dis","NuWro DIS", "f");
   legendb->Draw();
}

if (t>4)
{TLegend *legendb = new TLegend(0.65,0.6,0.9,0.75);
   legendb->AddEntry("t2k","T2K", "lep");
   legendb->AddEntry("ccqe","NuWro CCQE", "f");
   legendb->AddEntry("res","NuWro RES", "f");
   legendb->AddEntry("mec","NuWro MEC", "f");
   legendb->AddEntry("dis","NuWro DIS", "f");
   legendb->Draw();
}

  TPaveText *pt = new TPaveText(0.15,0.76,0.95,0.9, "brNDC");
  pt->AddText(napis4);
  pt->AddText(napis5);
  pt->AddText(napis6);
  //pt->AddText(napis7);
  //pt->AddText(napis8);
   pt->Draw();
  
   calphat->Print(napis2);
   delete hs, gan2, calphat, haqel, hares, hamec, hadis, pt;
   used2+=binning[t];
}


//////////////////////////
///////// stack (complete) histograms (mec)
///////////////////////////

pointer=0;
used2=0;
for (int t=0; t<cosbins; t++)//figures in cosine bins
{
  
for (int j=0; j<totalbins; j++)//identification of the first momentum bin
{
 if (aggr[j]==mombins*t)
   pointer=j;
}

  double a=cos_left[t];
  double b=cos_left[t+1];//start end of the cos bin, this part is easy
  int c = int (100*a+100);
  
  int mombins_new=binning[t]-1;//number of momentum bins
  double mom_left_new [mombins_new+1];//array with new binning
  
  mom_left_new[0]=5;//to avoid problems in log scale plots
  int mom_run=0;
  
  for (int s=1; s<mombins_new+1; s++)
  {mom_left_new[s]=mom_left[sklejka[pointer]+mom_run];
    mom_run+=sklejka[pointer+s];
  }
  
sprintf (napis, "sztuka_%d",t);
  sprintf (napis2, "T2K_CC0pi_ciro_anumu_%d_stack_complete_mec.png",c);
  //sprintf (napis3, "T2K CC0#pi CH, new, anumu, %.2f<cos#theta_{#mu}<%.2f, #chi^{2}=%.1f",a,b,chikwa[t]);
  sprintf (napis3, "T2K CC0#pi CH, new, anumu, %.2f<cos#theta_{#mu}<%.2f",a,b);
  sprintf (napis4, "kinematical characteristics of MEC events (in bins from left energy and momentum transfer)");
  
  sprintf (napis5, "1) w=%.0f#pm%.0f, ",w_average_mec[used2], w_sigma_mec[used2]);
  for (int s=1; s<mombins_new; s++)
sprintf(napis5 + strlen(napis5),"%d) w=%.0f#pm%.0f, ",s+1, w_average_mec[used2+s], w_sigma_mec[used2+s]);

 sprintf (napis6, "1) q=%.0f#pm%.0f, ",q_average_mec[used2], q_sigma_mec[used2]);
  for (int t=1; t<mombins_new; t++)
sprintf(napis6 + strlen(napis6),"%d) q=%.0f#pm%.0f, ",t+1, q_average_mec[used2+t], q_sigma_mec[used2+t]);

THStack *hs = new THStack("hs","Stacked 1D histograms");
hs->SetMinimum(0.0);

if (t==0)
hs->SetMaximum(4.9e-39);

if (t==1)
hs->SetMaximum(2e-39);

if (t==2)
hs->SetMaximum(6.5e-39);

if (t==3)
hs->SetMaximum(6.5e-39);

if (t==4)
hs->SetMaximum(8e-39);

if (t==5)
hs->SetMaximum(9e-39);

if (t==6)
hs->SetMaximum(10e-39);
 
if (t==7)
hs->SetMaximum(10e-39);
  
if (t==8)
hs->SetMaximum(6.25e-39);
   
TH1D *haqel = new TH1D ("haqel", "haqel", mombins_new, mom_left_new);
 TH1D *hares = new TH1D ("hares", "hares", mombins_new, mom_left_new);
 TH1D *hamec = new TH1D ("hamec", "hamec", mombins_new, mom_left_new);
 TH1D *hadis = new TH1D ("hadis", "hadis", mombins_new, mom_left_new);
 
double xdpt[mombins_new];
double ydpt[mombins_new];
double exdpt[mombins_new];
double eydpt[mombins_new];
 
for (int k=0; k<mombins_new; k++)
   {haqel->SetBinContent(k+1,cross_qel[pointer+k]);
     hares->SetBinContent(k+1,cross_res[pointer+k]);
     hamec->SetBinContent(k+1,cross_mec[pointer+k]);
     hadis->SetBinContent(k+1,cross_dis[pointer+k]);
     
    xdpt[k]=(mom_left_new[k+1]+mom_left_new[k])/2.0;
    exdpt[k]=(mom_left_new[k+1]-mom_left_new[k])/2.0;
    ydpt[k]= ciro_anu_sections[pointer+k];
    eydpt[k]=error[pointer+k];
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

TGraphErrors *gan2 = new TGraphErrors(mombins_new,xdpt,ydpt,exdpt,eydpt);
 

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
   calphat->SetLogx();

hs->Draw();
hs->SetTitle(napis3);
hs->SetName("phit");
hs->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
hs->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
hs->GetXaxis()->SetTitle("muon momentum [MeV]");
hs->GetYaxis()->SetTitle("cross section (per nucleon) [cm^{2}/GeV^{2}/nucleon]");
hs->Draw();

gan2->SetName("t2k");
gan2->SetLineWidth(3);
gan2->SetLineColor(kBlack);
gan2->Draw("E SAME"); 

if (t<5)
{TLegend *legendb = new TLegend(0.15,0.6,0.4,0.75);
   legendb->AddEntry("t2k","T2K", "lep");
   legendb->AddEntry("ccqe","NuWro CCQE", "f");
   legendb->AddEntry("res","NuWro RES", "f");
   legendb->AddEntry("mec","NuWro MEC", "f");
   legendb->AddEntry("dis","NuWro DIS", "f");
   legendb->Draw();
}

if (t>4)
{TLegend *legendb = new TLegend(0.65,0.6,0.9,0.75);
   legendb->AddEntry("t2k","T2K", "lep");
   legendb->AddEntry("ccqe","NuWro CCQE", "f");
   legendb->AddEntry("res","NuWro RES", "f");
   legendb->AddEntry("mec","NuWro MEC", "f");
   legendb->AddEntry("dis","NuWro DIS", "f");
   legendb->Draw();
}

  TPaveText *pt = new TPaveText(0.15,0.76,0.95,0.9, "brNDC");
  pt->AddText(napis4);
  pt->AddText(napis5);
  pt->AddText(napis6);
  //pt->AddText(napis7);
  //pt->AddText(napis8);
   pt->Draw();
  
   calphat->Print(napis2);
   delete hs, gan2, calphat, haqel, hares, hamec, hadis, pt;
   used2+=binning[t];
}

}

