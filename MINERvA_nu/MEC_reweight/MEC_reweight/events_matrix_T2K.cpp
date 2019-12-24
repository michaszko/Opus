#define Q_OMEGA_BINS 24
#define Q_OMEGA_RANGE 1200/Q_OMEGA_BINS

#define COS_BINS 9
#define MOM_BINS 14

double mec_events[Q_OMEGA_BINS][Q_OMEGA_BINS][58] = {0};
double cross_section = 0;
int good_events = 0;

double mom_left [MOM_BINS + 1] = {5, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 2000, 3000, 5000, 30000};
double cos_left [COS_BINS + 1] = { -1.0, 0.2, 0.6, 0.7, 0.8, 0.85, 0.9, 0.94, 0.98, 1.0};

bool get_cross_section(ifstream &file) {
    string ff;
    for (int i = 0; i < 9; ++i) {
        getline(file, ff);
    }
    file >> cross_section;
    if (int(cross_section) == 8) {
        file >> cross_section >> cross_section >> cross_section;
        cout << "Cross section of events: " << cross_section << endl;
    } else {
        return false;
    }

    return true;
}

int momentum_bin(event *e) {
    double muonmom = e->out[0].momentum();
    for (int mom_bin = 0; mom_bin < MOM_BINS; ++mom_bin) {
        if ((muonmom > mom_left[mom_bin]) && (muonmom <= mom_left[mom_bin + 1])) {
            return mom_bin;
        }
    }

    return -1;
}

int cosine_bin(event *e) {
    double muoncos = e->out[0].z / e->out[0].momentum();
    for (int cos_bin = 0; cos_bin < COS_BINS; ++cos_bin) {
        if ((muoncos > cos_left[cos_bin]) && (muoncos <= cos_left[cos_bin + 1])) {
            return cos_bin;
        }
    }

    return -1;
}


int final[COS_BINS][MOM_BINS] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 2, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    6, 7, 8, 9, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11,
    12, 13, 14, 15, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17,
    18, 19, 20, 21, 22, 22, 23, 23, 24, 24, 24, 24, 24, 24,
    25, 26, 27, 28, 29, 29, 30, 30, 31, 31, 32, 32, 32, 32,
    33, 33, 34, 35, 36, 36, 37, 37, 37, 38, 38, 39, 39, 39,
    40, 40, 41, 42, 43, 43, 44, 44, 45, 46, 47, 48, 49, 49,
    50, 50, 50, 51, 51, 52, 52, 53, 53, 54, 54, 55, 56, 57
};

void export_to_file() {
    ofstream python("Data/T2K_neutrino/events.py");

    python << "import numpy as np\n\n";
    python << "events = [\n";
    for (int q = 0; q < Q_OMEGA_BINS; ++q) {
        for (int omega = 0; omega < Q_OMEGA_BINS; ++omega) {
            python << "    ";
            for (int l = 0; l < 58; ++l) {
                python << mec_events[q][omega][l] << ", ";
            }
            python << "\n";
        }
    }

    python << "]\n\n";
    python << "events = np.array(events).reshape(" << Q_OMEGA_BINS << ", " << Q_OMEGA_BINS << ", " << 58 << ")\n\n";

    python << "nof_events = " << good_events << '\n';
    python << "cross_section = " << cross_section << '\n';

    python.close();
}

void events_matrix_T2K() {
    ifstream cross_file("Data/T2K_neutrino/neutrino_T2K_mec.txt");
    if (!get_cross_section(cross_file)) {
        cout << "Cannot retrieve cross section info, aborting." << endl;
        return;
    }

    TFile *data = new TFile("Data/T2K_neutrino/neutrino_T2K_mec.root");
    TTree *tree = (TTree *)data->Get("treeout");
    event *e = new event();
    tree->SetBranchAddress("e", &e);

    int events = tree->GetEntries();

    cout << "Analyzing " << events << " events" << endl;

    for (int k = 0; k < events; k++) {
        // break;
        cout << k  * 100.0 / events << "%\r";
        tree->GetEntry(k);

        if (!e->fof(211, -211, 111)) {
            int mom_bin = momentum_bin(e);
            int cos_bin = cosine_bin(e);

            double energy_transfer = e->in[0].t - e->out[0].t;
            double momentum_transfer = sqrt(pow(energy_transfer, 2) - e->q2());

            int q = energy_transfer / Q_OMEGA_RANGE;
            int omega = momentum_transfer / Q_OMEGA_RANGE;

            if (mom_bin >= 0 && cos_bin >= 0 &&  q < Q_OMEGA_BINS && omega < Q_OMEGA_BINS) {
                mec_events[q][omega][final[cos_bin][mom_bin]]++;
                good_events += 1;
            }
        }
    }

    int sklejka[58] = {
        14,
        1, 1, 1, 1, 10,
        1, 1, 1, 1, 2, 8,
        1, 1, 1, 1, 2, 8,
        1, 1, 1, 1, 2, 2, 6,
        1, 1, 1, 1, 2, 2, 2, 4,
        2, 1, 1, 2, 3, 2, 3,
        2, 1, 1, 2, 2, 1, 1, 1, 1, 2,
        3, 2, 2, 2, 2, 1, 1, 1
    };

    int suma[58];
    suma[0] = 0;
    for (int i = 1; i < 58; ++i) {
        suma[i] = suma[i - 1] + sklejka[i - 1];
    }

    for (int b = 0; b < 58; b++) {
        double norm = cross_section * 1e3 / events;
        norm /= (cos_left[suma[b] / 14 + 1] - cos_left[suma[b] / 14]);
        norm /= (mom_left[suma[b] % 14 + sklejka[b]] - mom_left[suma[b] % 14]);

        for (int q = 0; q < Q_OMEGA_BINS; ++q) {
            for (int w = 0; w < Q_OMEGA_BINS; ++w) {
                mec_events[q][w][b] *= norm;
            }
        }
    }


    export_to_file();

    printf("%d events exported\n", good_events);
}