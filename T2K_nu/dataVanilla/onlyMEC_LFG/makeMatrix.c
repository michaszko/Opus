const int mombins = 14;
const int cosbins = 9;
const int dim = 24;

int noPionsEventsCounter, x, y, xx, yy, row, column, events,binsBefore;
int matrix[600][600];

double MyPi = 3.14159265358979323846;

double mom_left[mombins + 1] = {5,    300,  400,  500,  600,  700,  800,  900, 1000, 1250, 1500, 2000, 3000, 5000, 30000};
double cos_left[cosbins + 1] = {-1,   0.2, 0.6,  0.7,  0.8, 0.85, 0.9, 0.94, 0.98, 1};

int binning[cosbins] = {1, 5, 6, 6, 7, 8, 7, 10, 8};

char fileName[100] = "onlyMEC.root";

ofstream mat("../../data/matrix_LFG.dat");

/////////////////////////////////////////////////////////////////////////

int position_cos(double a, double *taba, int tabaSize) {
  int xa = 0;

  for (int i = 0; i < tabaSize; i++) {
    if (a <= taba[i])
      break;

    xa = i;
  }
  return xa;
}

int position_mom(double a, double *taba, int tabaSize, int angle) {
  int xa = -1;

  for (int i = 0; i < binning[angle]; i++) {
    if (a <= taba[i])
      break;

    xa = i;
  }
  return xa;
}

///////////////////////////////////////

void makeMatrix() {
  TFile *Input = new TFile(fileName);
  TTree *tt1 = (TTree *)Input->Get("treeout");
  event *e = new event();
  tt1->SetBranchAddress("e", &e);

  events = tt1->GetEntries();

  for (int k = 0; k < events; k++) {
    cout << "\r"
         << "Events already\t" << k << "/" << events;

    tt1->GetEntry(k);

    if (e->fof(211, -211, 111) == 0 && e->fof(321, -321, 311) == 0 &&
        e->fof(-311, 130, 310) == 0) { // no pions
      noPionsEventsCounter++;

      double muoncos = e->out[0].z / e->out[0].momentum();
      double muonmom = e->out[0].momentum();
      double entrans = e->in[0].t - e->out[0].t;
      double momtrans2 = entrans * entrans - e->q2();
      double momtrans = sqrt(momtrans2);

      if (muonmom >= 5.00 && muonmom <= 30000.00) {
        x = position_cos(muoncos, cos_left, cosbins);
        y = position_mom(muonmom, mom_left, mombins, x);

        binsBefore = 0;

        for(int i = 0; i < x; i++)
        	binsBefore+=binning[i];

        xx = int(entrans / 1200.0 * dim);
        yy = int(sqrt(momtrans2) / 1200.0 * dim);

        // row = y + 13 * x ;
        row = binsBefore + y;
        column = xx * dim + yy;

        if(y != -1)
       		matrix[row][column]++;
      }
    }
  }

  cout << noPionsEventsCounter << "\n"
       << float(noPionsEventsCounter / events) << "\n";

  cout << "Now priting matrix to file...\n";

  for (int i = 0; i < 58; ++i) {
    for (int j = 0; j < dim * dim; ++j) {
      mat << matrix[i][j] << "\t";
    }
    mat << endl;
  }

  cout << "Finished!\n";
}