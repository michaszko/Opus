const int transbins_nu = 13;
const int longbins_nu = 12;

int noPionsEventsCounter = 0, events, inHistogram = 0;

double MyPi = 3.14159265358979323846;
double kosmin = cos(20.0 * MyPi / 180.0);

char fileNameTXT[100]="chi2_test.root.txt";
char fileNameROOT[100]="chi2_test.root";

ofstream cros_nuwro_total("../../data/cros_total_nuwro_test.dat");
ofstream cros_nuwro_mec("../../data/cros_mec_nuwro_test.dat");

double trans_nu_left[transbins_nu + 1] = {0, 75, 150, 250, 325, 400, 475, 550, 700, 850, 1000, 1250, 1500, 2500};
double long_nu_left[longbins_nu + 1] = {1500,  2000,  2500, 3000, 3500, 4000,  4500,  5000, 6000, 8000, 10000, 15000, 20000};

double results_nu[longbins_nu + 2][transbins_nu + 2];
double cross_results_nu[longbins_nu + 1][transbins_nu + 1];

double results_nu_qel[longbins_nu + 2][transbins_nu + 2];
double cross_results_nu_qel[longbins_nu + 1][transbins_nu + 1];

double results_nu_res[longbins_nu + 2][transbins_nu + 2];
double cross_results_nu_res[longbins_nu + 1][transbins_nu + 1];

double results_nu_mec[longbins_nu + 2][transbins_nu + 2];
double cross_results_nu_mec[longbins_nu + 1][transbins_nu + 1];

double results_nu_dis[longbins_nu + 2][transbins_nu + 2];
double cross_results_nu_dis[longbins_nu + 1][transbins_nu + 1];

/////////////////////////////////////

int position(double a, double *taba, int tabaSize) {
  int xa = 0;

  for (int i = 0; i < tabaSize; i++) {
    if (a <= taba[i]) break;

    xa = i;
  }
  return xa;
}

/////////////////////////////////////

void asia_a_daniel(TFile *Input,
                   double table_nu[longbins_nu + 2][transbins_nu + 2],
                   double table_nu_qel[longbins_nu + 2][transbins_nu + 2],
                   double table_nu_res[longbins_nu + 2][transbins_nu + 2],
                   double table_nu_mec[longbins_nu + 2][transbins_nu + 2],
                   double table_nu_dis[longbins_nu + 2][transbins_nu + 2]) {
  TH2D *h2 = new TH2D("h2", "h2", longbins_nu, long_nu_left, transbins_nu,trans_nu_left);
  TH2D *h2qel = new TH2D("h2qel", "h2qel", longbins_nu, long_nu_left,transbins_nu, trans_nu_left);
  TH2D *h2res = new TH2D("h2res", "h2res", longbins_nu, long_nu_left,transbins_nu, trans_nu_left);
  TH2D *h2mec = new TH2D("h2mec", "h2mec", longbins_nu, long_nu_left,transbins_nu, trans_nu_left);
  TH2D *h2dis = new TH2D("h2dis", "h2dis", longbins_nu, long_nu_left,transbins_nu, trans_nu_left);

  TTree *tt1 = (TTree *)Input->Get("treeout");
  event *e = new event();
  tt1->SetBranchAddress("e", &e);

  events = tt1->GetEntries();

  for (int k = 0; k < events; k++) {
    cout << "\r" << "Events already\t" << k << "/" << events;

    tt1->GetEntry(k);

    if (e->fof(211, -211, 111) == 0 && e->fof(321, -321, 311) == 0 && e->fof(-311, 130, 310) == 0)  // no pions
    {
      double kos = e->out[0].z / e->out[0].momentum();

      if (kos > kosmin) {
      	noPionsEventsCounter++;

        double muonlong = e->out[0].z;
        double muontrans = sqrt((e->out[0].x) * (e->out[0].x) + (e->out[0].y) * (e->out[0].y));
        double entrans = e->in[0].t - e->out[0].t;
        double momtrans2 = entrans * entrans - e->q2();

        if (entrans > 1200 or sqrt(momtrans2) > 1200) continue;

        if (muonlong >= 1500.00 && muonlong <= 20000.00) inHistogram++;

	        h2->Fill(muonlong, muontrans);

	        if (e->flag.qel == 1) h2qel->Fill(muonlong, muontrans);

	        if (e->flag.res == 1) h2res->Fill(muonlong, muontrans);

	        if (e->flag.mec == 1) h2mec->Fill(muonlong, muontrans);

	        if (e->flag.dis == 1) h2dis->Fill(muonlong, muontrans);
      }
    }
  }

  cout << endl << "No pions events and cosinus (in percent)" << (double)noPionsEventsCounter / (double)events * 100 << "%" << endl;

  cout << endl << "In in histogram " << inHistogram << endl;

  delete e;

  for (int j = 0; j < longbins_nu + 2; j++) {
    for (int k = 0; k < transbins_nu + 2; k++) {
      table_nu[j][k]      = h2->GetBinContent(j, k);
      table_nu_qel[j][k]  = h2qel->GetBinContent(j, k);
      table_nu_res[j][k]  = h2res->GetBinContent(j, k);
      table_nu_mec[j][k]  = h2mec->GetBinContent(j, k);
      table_nu_dis[j][k]  = h2dis->GetBinContent(j, k);
    }
  }
  delete h2, h2qel, h2res, h2mec, h2dis;
}

/////////////////////////////////////

double norm_cc(ifstream &Input) {
  string ff;
  double answer = 0;
  double aa;

  while (!Input.eof()) {
    getline(Input, ff);

    for (int kk = 0; kk < 10; kk++) {
      for (int jj = 0; jj < 4; jj++) {
        Input >> aa;
      }

      if (kk == 0) answer += aa;
      if (kk == 2) answer += aa;
      if (kk == 4) answer += aa;
      if (kk == 6) answer += aa;
      if (kk == 8) answer += aa;
    }
  }
  return answer;
}

//////////////////////////////////////////
/////////////////////////////////////////
void dataProcessing()
/////////////////////////////////////////
////////////////////////////////////////
{
  cout << "Loading file " << endl;

  ifstream Input2 (fileNameTXT);
  if (!Input2.is_open())
    cout << "Huston, we have a problem with loading file :| " << endl;
  else
    cout << "File .txt loaded succesfuly" << endl;

  cout << "Norming now" << endl;
  double cross_minerva_nu = norm_cc(Input2);
  cout << "lol" << endl;

  TFile *tf2 = new TFile(fileNameROOT);
  if (!tf2->IsOpen())
    cout << "Huston, we have a problem with loading file :| " << endl;
  else
    cout << "File .root loaded succesfuly" << endl;

  cout << "asia_a_daniel" << endl;

  asia_a_daniel(tf2, results_nu, results_nu_qel, results_nu_res, results_nu_mec, results_nu_dis);

  cout << "\nFilling files" << endl;

  for (int j = 1; j < longbins_nu + 1; j++) {
    for (int k = 1; k < transbins_nu + 1; k++) {
      cross_results_nu[j][k] = 
        results_nu[j][k] / events * cross_minerva_nu /
        (long_nu_left[j] - long_nu_left[j - 1]) /
        (trans_nu_left[k] - trans_nu_left[k - 1]) * 1e6 * 12 / 13;

      cross_results_nu_qel[j][k] =
          results_nu_qel[j][k] / events * cross_minerva_nu /
          (long_nu_left[j] - long_nu_left[j - 1]) /
          (trans_nu_left[k] - trans_nu_left[k - 1]) * 1e6 * 12 / 13;

      cross_results_nu_res[j][k] =
          results_nu_res[j][k] / events * cross_minerva_nu /
          (long_nu_left[j] - long_nu_left[j - 1]) /
          (trans_nu_left[k] - trans_nu_left[k - 1]) * 1e6 * 12 / 13;

      cross_results_nu_mec[j][k] =
          results_nu_mec[j][k] / events * cross_minerva_nu /
          (long_nu_left[j] - long_nu_left[j - 1]) /
          (trans_nu_left[k] - trans_nu_left[k - 1]) * 1e6 * 12 / 13;

      cross_results_nu_dis[j][k] =
          results_nu_dis[j][k] / events * cross_minerva_nu /
          (long_nu_left[j] - long_nu_left[j - 1]) /
          (trans_nu_left[k] - trans_nu_left[k - 1]) * 1e6 * 12 / 13;

      cros_nuwro_total  << cross_results_nu[j][k] << "\t";
      cros_nuwro_mec    << cross_results_nu_mec[j][k] << "\t";
    }
    cros_nuwro_total  << endl;
    cros_nuwro_mec    << endl;
  }
}
