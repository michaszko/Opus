// #include "ciro_cov.h"
// #include "ciro_data_new.h"
// #include "mec_scaling_april23_corr.h"
// #include "mec_scaling_t2k_nu_april23_corr.h"

const int mombins = 14;
const int cosbins = 9;
double mom_left[mombins + 1] = { 5, 300, 400, 500, 600, 700, 800, 900,
    1000, 1250, 1500, 2000, 3000, 5000, 30000 };
double cos_left[cosbins + 1] = { -1, 0.2, 0.6, 0.7, 0.8,
    0.85, 0.9, 0.94, 0.98, 1 };

double results[mombins + 2][cosbins + 2];
double results_qel[mombins + 2][cosbins + 2];
double results_res[mombins + 2][cosbins + 2];
double results_mec[mombins + 2][cosbins + 2];
double results_dis[mombins + 2][cosbins + 2];

const int totalbins = 58;
double cross[totalbins];
double cross_qel[totalbins];
double cross_res[totalbins];
double cross_mec[totalbins];
double cross_dis[totalbins];

int binning[cosbins] = { 1, 5, 6, 6, 7, 8, 7, 10, 8 };

int sklejka[totalbins] = { 14, 1, 1, 1, 1, 10, 1, 1, 1, 1, 2, 8, 1, 1, 1,
    1, 2, 8, 1, 1, 1, 1, 2, 2, 6, 1, 1, 1, 1, 2,
    2, 2, 4, 2, 1, 1, 2, 3, 2, 3, 2, 1, 1, 2, 2,
    1, 1, 1, 1, 2, 3, 2, 2, 2, 2, 1, 1, 1 };

///////////////////////////////////////////////////////
double mscaling(double momtran, double entran)
{
    int momnumber = int(momtran / 100);
    int ennumber = int(entran / 100);

    // cout<<momtran<<"  "<<momnumber<<"  "<<entran<<"  "<<ennumber<<"
    // "<<mec_scaling [momnumber][ennumber]<<endl;  return mec_scaling
    // [momnumber][ennumber];  return mec_scaling_t2k_nu [momnumber][ennumber];
    return 1;
}
//////////////////////////////////////////////////////

ofstream nuwro("1709_t2k_CH_ciro_numu_report.dat");
ofstream nuwro2("cross_section.txt");

void asia_a(TFile* Input, int ilezdarzen,
    double table[mombins + 2][cosbins + 2],
    double table_qel[mombins + 2][cosbins + 2],
    double table_res[mombins + 2][cosbins + 2],
    double table_mec[mombins + 2][cosbins + 2],
    double table_dis[mombins + 2][cosbins + 2],)
{
    TH2D* h1 = new TH2D("h1", "h1", mombins, mom_left, cosbins, cos_left);
    TH2D* h1qel = new TH2D("h1qel", "h1qel", mombins, mom_left, cosbins, cos_left);
    TH2D* h1res = new TH2D("h1res", "h1res", mombins, mom_left, cosbins, cos_left);
    TH2D* h1mec = new TH2D("h1mec", "h1mec", mombins, mom_left, cosbins, cos_left);
    TH2D* h1dis = new TH2D("h1dis", "h1dis", mombins, mom_left, cosbins, cos_left);

    TTree* tt1 = (TTree*)Input->Get("treeout");
    event* e = new event();
    tt1->SetBranchAddress("e", &e);

    events = tt1->GetEntries();

    for (int k = 0; k < events; k++) {
        tt1->GetEntry(k);
        {
            if (e->fof(211, -211, 111) == 0) // no pions
            {
                double muoncos = e->out[0].z / e->out[0].momentum();
                double muonmom = e->out[0].momentum();
                double entrans = e->in[0].t - e->out[0].t;
                double momtrans2 = entrans * entrans - e->q2();
                double momtrans = sqrt(momtrans2);

                if (e->flag.mec == 1)
                    h1->Fill(muonmom, muoncos, mscaling(momtrans, entrans));
                else
                    h1->Fill(muonmom, muoncos);

                if (e->flag.qel == 1) {
                    h1qel->Fill(muonmom, muoncos);
                }

                if (e->flag.res == 1)
                    h1res->Fill(muonmom, muoncos);

                if (e->flag.mec == 1) {
                    h1mec->Fill(muonmom, muoncos, mscaling(momtrans, entrans));
                }

                if (e->flag.dis == 1)
                    h1dis->Fill(muonmom, muoncos);
            }
        }
    }
    // delete e;

    for (int j = 0; j < mombins + 2; j++) {
        for (int k = 0; k < cosbins + 2; k++) {
            table[j][k] = h1->GetBinContent(j, k);
            table_qel[j][k] = h1qel->GetBinContent(j, k);
            table_res[j][k] = h1res->GetBinContent(j, k);
            table_mec[j][k] = h1mec->GetBinContent(j, k);
            table_dis[j][k] = h1dis->GetBinContent(j, k);
        }
    }
    // delete h1, h1qel, h1res, h1mec, h1dis, h1qel_w, h1qel_w2, h1qel_q, h1qel_q2,
    // h1mec_w, h1mec_w2, h1mec_q, h1mec_q2;
}

double norm_cc(ifstream& Input)
{
    string ff;
    double answer = 0;
    double aa;

    while (!Input.eof()) {
        getline(Input, ff);
        for (int kk = 0; kk < 10; kk++) {
            for (int jj = 0; jj < 4; jj++) {
                Input >> aa;
            }
            if (kk == 0)
                answer += aa;
            if (kk == 2)
                answer += aa;
            if (kk == 4)
                answer += aa;
            if (kk == 6)
                answer += aa;
            if (kk == 8)
                answer += aa;
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
    ifstream Input("a_t2k_cc_CH_lfg_500kilo.root.txt");
    double cross_t2k = norm_cc(Input);

    TFile* tf1 = new TFile("a_t2k_cc_CH_lfg_500kilo.root");
    asia_a(tf1, events, results, results_qel, results_res, results_mec, results_dis);
    delete tf1;

    int aggr[totalbins]; // how many bins before

    aggr[0] = 0;
    for (int j = 1; j < totalbins; j++)
        aggr[j] = aggr[j - 1] + sklejka[j - 1];

    ////////////////////////////////
    //////// cross sections, errors
    ////////////////////////////////

    nuwro << "cross section results" << endl;
    for (int j = 0; j < totalbins; j++) {
        int dlugosc = sklejka[j];
        int koss = aggr[j] / mombins; // which cos bin
        int momm = aggr[j] - koss * mombins; // which mom bin
        double bin_start = mom_left[momm];
        double bin_end = mom_left[momm + dlugosc]; // correct understanding of binning

        double sumka = 0;
        double sumka_qel = 0;
        double sumka_res = 0;
        double sumka_mec = 0;
        double sumka_dis = 0;

        for (int s = momm + 1; s < momm + dlugosc + 1; s++) {
            sumka += results[s][koss + 1];
            sumka_qel += results_qel[s][koss + 1];
            sumka_res += results_res[s][koss + 1];
            sumka_mec += results_mec[s][koss + 1];
            sumka_dis += results_dis[s][koss + 1];
        }

        cross[j] = sumka / events * cross_t2k / (cos_left[koss + 1] - cos_left[koss]) / (mom_left[momm + dlugosc] - mom_left[momm]) * 1e3;
        cross_qel[j] = sumka_qel / events * cross_t2k / (cos_left[koss + 1] - cos_left[koss]) / (mom_left[momm + dlugosc] - mom_left[momm]) * 1e3;
        cross_res[j] = sumka_res / events * cross_t2k / (cos_left[koss + 1] - cos_left[koss]) / (mom_left[momm + dlugosc] - mom_left[momm]) * 1e3;
        cross_mec[j] = sumka_mec / events * cross_t2k / (cos_left[koss + 1] - cos_left[koss]) / (mom_left[momm + dlugosc] - mom_left[momm]) * 1e3;
        cross_dis[j] = sumka_dis / events * cross_t2k / (cos_left[koss + 1] - cos_left[koss]) / (mom_left[momm + dlugosc] - mom_left[momm]) * 1e3;
        // nuwro<<cross[j]<<"  "<<cross_qel[j]<<"  "<<cross_res[j]<<" "<<cross_mec[j]<<"  "<<cross_dis[j]<<"  ";
        nuwro2 << cross[j] << " ";
    }

    double error[totalbins];
    nuwro << "cross sections and errors" << endl;

    for (int s = 0; s < totalbins; s++) {
        error[s] = // sqrt( //ciro_nu_cov_norm [s*(totalbins+1)]+ ciro_nu_cov_shape
            // [s*(totalbins+1)] +
            ciro_nu_errors[s];
        nuwro << cross[s] << "  " << ciro_nu_sections[s] << "  " << error[s]
              << endl;
    }
    nuwro << endl;
}
