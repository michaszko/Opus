const int transbins_nu = 13;
const int longbins_nu = 12;
const int dim = 24;

int noPionsEventsCounter, x, y, xx, yy, row, column, events;
int matrix[600][600];

double MyPi = 3.14159265358979323846;
double kosmin = cos(20.0 * MyPi / 180.0);

double trans_nu_left[transbins_nu + 1] = {0, 75, 150, 250, 325, 400, 475, 550, 700, 850, 1000, 1250, 1500, 2500};
double long_nu_left[longbins_nu + 1] = {1500,  2000,  2500, 3000, 3500, 4000, 4500,  5000, 6000, 8000, 10000, 15000, 20000};

char fileName[100]="chi2_test.root";

ofstream mat("../../data/matrix_test.dat");

/////////////////////////////////////////////////////////////////////////

int position(double a, double *taba, int tabaSize) {
  int xa = 0;

  for (int i = 0; i < tabaSize; i++) {
    if (a <= taba[i]) break;

    xa = i;
  }
  return xa;
}

///////////////////////////////////////

void makeMatrix()
{
	TFile *Input = new TFile(fileName);
	TTree *tt1 = (TTree *)Input->Get("treeout");
	event *e = new event();
	tt1->SetBranchAddress("e", &e);

	events =  tt1->GetEntries();

	for (int k = 0; k < events; k++){
		cout << "\r" << "Events already\t" << k << "/" << events;

		tt1->GetEntry(k);

		if (e->flag.mec == 1){

			if (e->fof(211, -211, 111) == 0 && e->fof(321, -321, 311) == 0 && e->fof(-311, 130, 310) == 0){ // no pions
			  double kos = e->out[0].z / e->out[0].momentum();

			  if (kos > kosmin) {
			  	noPionsEventsCounter++;

			    double muonlong = e->out[0].z;
			    double muontrans = sqrt((e->out[0].x) * (e->out[0].x) + (e->out[0].y) * (e->out[0].y));
			    double entrans = e->in[0].t - e->out[0].t;
			    double momtrans2 = entrans * entrans - e->q2();

			    if (muonlong >= 1500.00 && muonlong <= 20000.00) {
					x = position(muonlong, long_nu_left, longbins_nu);
					y = position(muontrans, trans_nu_left, transbins_nu);

		      		xx = int(entrans/1200.0 * dim);
					yy = int(sqrt(momtrans2)/1200.0 * dim);

					row = y + 13 * x ;
					column = xx * dim + yy;

					matrix[row][column]++;
			    }
			  }
			}
		}
	}

	cout << endl << noPionsEventsCounter << "\n" << (float)noPionsEventsCounter/(float)events * 100 << "%" << "\n";

	cout << "Now priting matrix to file...\n";

	for (int i = 0; i < 156; ++i){
		for (int j = 0; j < dim*dim; ++j){
				mat << matrix[i][j] << "\t";
			}
			mat << endl;
	}

	cout << "Finished!\n";
}