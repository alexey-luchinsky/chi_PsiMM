#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVector3R.hh"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TF1.h"
using namespace std;

void write_histogram_to_file(TH1F &histogram, string file_name) {
    const char *__file_name__ = file_name.c_str();
    remove(__file_name__);
    ofstream file;
    cout << "Extporting histogram " << file_name << endl;
    file.open(__file_name__);
    for (int i = 1; i <= histogram.GetNbinsX(); i++)
        file << setiosflags(ios::scientific) << histogram.GetBinCenter(i) <<
        " " << setiosflags(ios::scientific) << histogram.GetBinContent(i) / histogram.GetBinWidth(i) <<
        " " << setiosflags(ios::scientific) << histogram.GetBinError(i) / histogram.GetBinWidth(i) << endl;

    file.close();
}

void TLorentzVector_2_EvtRector4R(EvtVector4R &out, TLorentzVector *in) {
    out.set(in->T(), in->X(), in->Y(), in->Z());
}

double getPolarCos(EvtVector4R pPsi, EvtVector4R kk1);

int main(int argc, char const *argv[]) {
    const char *root_file_name;
    if (argc > 1)
        root_file_name = argv[1];
    else {
        cout << "format: ./read_data.exe inFileName [postFix=''] [nBins=20]" << endl;
        return 1;
    };
    string postFix="";
    if(argc>2) postFix=string(argv[2]);
    int nBins=20;
    if(argc>3) nBins=atoi(argv[3]);
    
    TFile file(root_file_name);
    TTree *moms = (TTree*) file.Get("moms");
    cout << moms->GetEntries() << " entries in the tree" << endl;
    TLorentzVector *_pPsi, *_k1, *_k2, *_kk1, *_kk2;
    double q2;
    moms->SetBranchAddress("q2", &q2);
    moms->SetBranchAddress("pPsi", &_pPsi);
    moms->SetBranchAddress("k1", &_k1);
    moms->SetBranchAddress("k2", &_k2);
    moms->SetBranchAddress("kk1", &_kk1);
    moms->SetBranchAddress("kk2", &_kk2);

    double minQ2 = moms->GetMinimum("q2"), maxQ2 = moms->GetMaximum("q2");
    //    int nBins = 20;
    //    TH1F *hQ2 = new TH1F("hQ2", "hQ2", nBins, minQ2, maxQ2);
    //    for (int iEv = 0; iEv < moms->GetEntries(); ++iEv) {
    //        moms->GetEntry(iEv);
    //        hQ2->Fill(q2);
    //    };
    //
    //    write_histogram_to_file(*hQ2, "hh.hst");


    // polarizations
    ofstream polFile;
    polFile.open(("pol"+postFix+".hst").c_str());
    EvtVector4R pPsi, k1, k2, kk1, kk2;
    TH1F *hPcos = new TH1F("hPcos", "hPcos", 20, -1, 1);
    hPcos->Sumw2();
    double dQ2 = (maxQ2 - minQ2) / nBins;
    for (double Q20 = minQ2 + dQ2 / 2; Q20 <= maxQ2; Q20 += dQ2) {
        hPcos->Clear();
        for (int iEv = 0; iEv < moms->GetEntries(); ++iEv) {
            moms->GetEntry(iEv);
            if (q2 < Q20 - dQ2 / 2 || q2 > Q20 + dQ2 / 2) continue;
            TLorentzVector_2_EvtRector4R(pPsi, _pPsi);
            TLorentzVector_2_EvtRector4R(kk1, _kk1);
            TLorentzVector_2_EvtRector4R(kk2, _kk2);
            kk1.applyBoostTo(pPsi, true);
            kk2.applyBoostTo(pPsi, true);
            double polarCos = getPolarCos(pPsi, kk1);
            hPcos->Fill(polarCos);
            if (iEv < 0 || abs(polarCos) > 1) {
                cout << " pPsi=" << pPsi << "  " << pPsi.mass() << endl;
                cout << " kk1=" << kk1 << "  " << kk1.mass() << endl;
                cout << " kk2=" << kk2 << "  " << kk2.mass() << endl;
                cout << " polarCos=" << getPolarCos(pPsi, kk1) << endl;
                cout << "========" << endl;
            };
        };
        //    write_histogram_to_file(*hPcos, "hPcos.hst");
        // fit it
        TF1 *myFit = new TF1("myFit", "[0]*(1+[1]*x*x)", 0, 1);
        myFit->SetParameter(0, 1);
        myFit->SetParameter(1, 0);
        hPcos->Fit("myFit", "q");
        TF1 *fit = hPcos->GetFunction("myFit");
        polFile << Q20 << " " << fit->GetParameter(1) << " " << fit->GetParError(1) << endl;
        cout << Q20 << " " << fit->GetParameter(1) << " " << fit->GetParError(1) << endl;
    };
    polFile.close();
    return 0;
}

double getPolarCos(EvtVector4R pPsi, EvtVector4R kk1) {
    EvtVector3R nPsi(pPsi.get(1), pPsi.get(2), pPsi.get(3));
    EvtVector3R nkk(kk1.get(1), kk1.get(2), kk1.get(3));
    return (nkk * nPsi) / nPsi.d3mag() / nkk.d3mag();
};
