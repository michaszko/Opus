#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std; 

double xSec(const char* fileName, bool all=false)
{
  int dyn, n;
  double ratio, sigma ,total=0;
  string comment;
 
  ifstream Input (fileName);

  getline (Input,comment); // line with headers

  while(Input>> dyn >> n>> ratio >>sigma) // lines with numbers
  {  
      if( all || dyn==0 || dyn==2 || dyn==4 || dyn==6 || dyn==8) 
        total += sigma;
  }
  
  return total;

}

char first[130];
char second[130];
char third[130];
char fourth[130];

ofstream file ("paris_2019.dat");

const double bins=20;

double en[20];
double cross_numu0[20];
double cross_numu1[20];
double cross_anumu0[20];
double cross_anumu1[20];

void mecscaling_sections()
{
  
  //cout<<xSec("nieves0numu_250.root.txt")<<endl;
  
  int ile=0;
  for (int E=250; E<5001; E+=250)
  {
	en[ile]=E;

	sprintf(first, "nieves0numu_%d.root.txt",E);
    	//ifstream Input1 (first);
	cross_numu0[ile]=xSec(first);
	//cout<<xSec(first)<<endl;

	sprintf(second, "nieves1numu_%d.root.txt",E);
    	//ifstream Input2 (second);
	cross_numu1[ile]=xSec(second);

	sprintf(third, "nieves0anumu_%d.root.txt",E);
    	//ifstream Input3 (third);
	cross_anumu0[ile]=xSec(third);

	sprintf(fourth, "nieves1anumu_%d.root.txt",E);
    	//ifstream Input4 (fourth);
	cross_anumu1[ile]=xSec(fourth);

    cout<<E<<"  "<<ile<<"  "<<cross_numu0[ile]<<"  "<<cross_numu1[ile]<<"  "<<cross_anumu0[ile]<<"  "<<cross_anumu1[ile]<<endl;
    ile++;
  }
  
TGraph *m1 = new TGraph (bins, en, cross_numu0);
TGraph *m2 = new TGraph (bins, en, cross_numu1);
TGraph *m3 = new TGraph (bins, en, cross_anumu0);
TGraph *m4 = new TGraph (bins, en, cross_anumu1);

////////////////////////
////////////////////////


TCanvas *c1cross = new TCanvas("c1cross","Test2",200,10,700,500);

   c1cross->Range(-0.2509977,-0.1717173,1.254176,1.076966);
   c1cross->SetFillColor(0);
   c1cross->SetBorderSize(2);
   c1cross->SetTickx();
   c1cross->SetTicky();
   c1cross->SetLeftMargin(0.1303828);
   c1cross->SetRightMargin(0.16937799);//originally 0.06
   c1cross->SetTopMargin(0.0612855);
   c1cross->SetBottomMargin(0.1375187);
   c1cross->SetFrameBorderMode(0);
   c1cross->SetFrameBorderMode(0);
   
   m1->SetTitle();
m1->SetName("m1");
m1->SetLineColor(3);
m1->SetLineWidth(5);
m1->SetLineStyle(3);

m2->SetName("m2");
m2->SetLineColor(2);
m2->SetLineWidth(5);
m2->SetLineStyle(2);

TMultiGraph *mg = new TMultiGraph();
mg->Add(m1);
mg->Add(m2);

mg->Draw("AC");

mg->GetXaxis()->SetTitle("Neutrino energy [MeV]");
mg->GetXaxis()->CenterTitle();
mg->GetYaxis()->SetTitle("cross section (per nucleon)");
mg->GetYaxis()->CenterTitle();

TLegend *leg = new TLegend(0.4,0.2,0.82,0.4);
   //leg->SetHeader("The Legend Title");
   //leg->AddEntry(h1,"Histogram filled with random numbers","f");
   leg->SetTextSize(0.04);
   leg->AddEntry("m1","Nieves model","l");
   leg->AddEntry("m2","Phenomenological model","l");
   leg->Draw();
  /* 
pt2 = new TPaveText(0.4, 0.15, 0.8, 0.27, "NDC"); // NDC sets coords
                                              // relative to pad dimensions
pt2->SetFillColor(0);
pt2->SetTextSize(0.08); 
text = pt2->AddText("PRELIMINARY");
pt2->Draw("SAME");       //to draw your text object
*/
TPaveText *pt3 = new TPaveText(0.2, 0.94, 0.8, 0.99, "NDC"); // NDC sets coords
                                              // relative to pad dimensions
pt3->SetFillColor(0);
pt3->SetTextSize(0.04); 
pt3->AddText("Neutrino");
pt3->Draw("SAME");       //to draw your text object

c1cross->Update();

c1cross->Print("mec_sections1_paris.png");

/*
 * m1->SetTitle();
m1->SetMarkerColor(1);
m1->GetXaxis()->SetTitle("Neutrino energy [MeV]");
m1->GetXaxis()->CenterTitle();
m1->GetYaxis()->SetTitle("cross section (per nucleon)");
m1->GetYaxis()->CenterTitle();
//m1->SetStats(kFALSE);
//m1->Draw("Box");
gStyle->SetPalette(1);
//hresult->SetContour(10);
m1->Draw("COLZ");
*/

/*
pt2 = new TPaveText(0.2, 0.94, 0.8, 0.99, "NDC"); // NDC sets coords
                                              // relative to pad dimensions
pt2->SetFillColor(0);
pt2->SetTextSize(0.04); 
text = pt2->AddText("MEC models cross sections (in NuWro)");
pt2->Draw("SAME");       //to draw your text object
c1cross->Update();
*/

//c1cross->Print("mec_sections.png");
////////////////////////
///////////////////////
TCanvas *c2cross = new TCanvas("c2cross","Test2",200,10,700,500);

   c2cross->Range(-0.2509977,-0.1717173,1.254176,1.076966);
   c2cross->SetFillColor(0);
   c2cross->SetBorderSize(2);
   c2cross->SetTickx();
   c2cross->SetTicky();
   c2cross->SetLeftMargin(0.1303828);
   c2cross->SetRightMargin(0.16937799);//originally 0.06
   c2cross->SetTopMargin(0.0612855);
   c2cross->SetBottomMargin(0.1375187);
   c2cross->SetFrameBorderMode(0);
   c2cross->SetFrameBorderMode(0);
   
   m3->SetTitle();
m3->SetName("m3");
m3->SetLineColor(3);
m3->SetLineWidth(5);
m3->SetLineStyle(3);

m4->SetName("m4");
m4->SetLineColor(2);
m4->SetLineWidth(5);
m4->SetLineStyle(2);

TMultiGraph *mg2 = new TMultiGraph();
mg2->Add(m3);
mg2->Add(m4);

mg2->Draw("AC");

mg2->GetXaxis()->SetTitle("Antineutrino energy [MeV]");
mg2->GetXaxis()->CenterTitle();
mg2->GetYaxis()->SetTitle("cross section (per nucleon)");
mg2->GetYaxis()->CenterTitle();

leg = new TLegend(0.4,0.2,0.82,0.4);
leg->SetTextSize(0.04);
   //leg->SetHeader("The Legend Title");
   //leg->AddEntry(h1,"Histogram filled with random numbers","f");
   leg->AddEntry("m3","Nieves model","l");
   leg->AddEntry("m4","Phenomenological model","l");
   leg->Draw();
   /*
pt4 = new TPaveText(0.4, 0.15, 0.8, 0.27, "NDC"); // NDC sets coords
                                              // relative to pad dimensions
pt4->SetFillColor(0);
pt4->SetTextSize(0.08); 
text = pt4->AddText("PRELIMINARY");
pt4->Draw("SAME");       //to draw your text object
*/
TPaveText *pt5 = new TPaveText(0.2, 0.94, 0.8, 0.99, "NDC"); // NDC sets coords
                                              // relative to pad dimensions
pt5->SetFillColor(0);
pt5->SetTextSize(0.04); 
pt5->AddText("Antineutrino");
pt5->Draw("SAME");       //to draw your text object


c2cross->Update();

c2cross->Print("mec_sections2_paris.png");


  
}
