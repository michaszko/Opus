#include "consts.h"

#define Q_OMEGA_BINS 24
#define Q_OMEGA_RANGE 1200/Q_OMEGA_BINS

#define TRANS_BINS 6
#define LONG_BINS 10

unsigned int mec_events[Q_OMEGA_BINS][Q_OMEGA_BINS][LONG_BINS][TRANS_BINS] = {0};
double cross_section = 0;
int good_events = 0;

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

int trans_momentum_bin(event *e) {
    double trans_momentum = sqrt( pow(e->out[0].x, 2) + pow(e->out[0].y, 2));
    for (int m_trans = 0; m_trans < TRANS_BINS; ++m_trans) {
        if ((trans_momentum > trans_nu_left[m_trans]) && (trans_momentum <= trans_nu_left[m_trans + 1])) {
            return m_trans;
        }
    }

    return -1;
}

int long_momentum_bin(event *e) {
    double long_momentum = e->out[0].z;
    for (int m_long = 0; m_long < LONG_BINS; ++m_long) {
        if ((long_momentum > long_nu_left[m_long]) && (long_momentum <= long_nu_left[m_long + 1])) {
            return m_long;
        }
    }

    return -1;
}

void export_to_file() {
    ofstream cpp("Data/antyneutrino_minerva/events.h");

    cpp << "int events[" << Q_OMEGA_BINS << "][" << Q_OMEGA_BINS << "][" << LONG_BINS << "][" << TRANS_BINS << "] = {\n";

    for (int q = 0; q < Q_OMEGA_BINS; ++q) {
        for (int omega = 0; omega < Q_OMEGA_BINS; ++omega) {
            for (int l = 0; l < LONG_BINS; ++l) {
                cpp << "    ";
                for (int t = 0; t < TRANS_BINS; ++t) {
                    cpp << mec_events[q][omega][l][t] << ", ";
                }
                cpp << "\n";
            }
        }
    }
    cpp << "};\n\n";

    cpp << "int nof_events = " << good_events << ";\n";
    cpp << "double cross_section = " << cross_section << ";\n";

    cpp.close();

    ofstream python("Data/antyneutrino_minerva/events.py");

    python << "import numpy as np\n\n";
    python << "events = [\n";
    for (int q = 0; q < Q_OMEGA_BINS; ++q) {
        for (int omega = 0; omega < Q_OMEGA_BINS; ++omega) {
            for (int l = 0; l < LONG_BINS; ++l) {
                python << "    ";
                for (int t = 0; t < TRANS_BINS; ++t) {
                    python << mec_events[q][omega][l][t] << ", ";
                }
                python << "\n";
            }
        }
    }

    python << "]\n\n";
    python << "events = np.array(events).reshape(" << Q_OMEGA_BINS << ", " << Q_OMEGA_BINS << ", " << LONG_BINS << ", " << TRANS_BINS << ")\n\n";

    python << "nof_events = " << good_events << '\n';
    python << "cross_section = " << cross_section << '\n';

    python.close();
}

void events_matrix() {
    ifstream cross_file("Data/antyneutrino_minerva/antyneutrino_minerva_mec.txt");
    if (!get_cross_section(cross_file)) {
        cout << "Cannot retrieve cross section info, aborting." << endl;
        return;
    }

    TFile *data = new TFile("Data/antyneutrino_minerva/antyneutrino_minerva_mec.root");
    TTree *tree = (TTree *)data->Get("treeout");
    event *e = new event();
    tree->SetBranchAddress("e", &e);

    int events = tree->GetEntries();

    cout << "Analyzing " << events << " events" << endl;

    for (int k = 0; k < events; k++) {
        tree->GetEntry(k);

        if (!e->fof(211, -211, 111) &&
                !e->fof(321, -321, 311) &&
                !e->fof(-311, 130, 310) &&
                e->out[0].z / e->out[0].momentum() > kosmin) {


            int trans_bin = trans_momentum_bin(e);
            int long_bin = long_momentum_bin(e);

            double energy_transfer = e->in[0].t - e->out[0].t;
            double momentum_transfer = sqrt(pow(energy_transfer, 2) - e->q2());

            int q = energy_transfer / Q_OMEGA_RANGE;
            int omega = momentum_transfer / Q_OMEGA_RANGE;

            // printf("%f, %f, %f, %f\n", trans_bin, long_bin, q, omega);

            if (long_bin >= 0 && trans_bin >= 0 &&  q < Q_OMEGA_BINS && omega < Q_OMEGA_BINS) {
                mec_events[q][omega][long_bin][trans_bin]++;
                good_events += 1;
            }
        }
    }

    export_to_file();

    printf("%d events exported\n", good_events);
}