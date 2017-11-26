#include <iostream>
#include <iomanip>
#include <fstream>

#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"

#include "kinematics.h"
#include "matr2.h"

using namespace std;

extern "C" {
    void init_(void);
    double ram2_(double &ecm, double *xm, double *p1, double *p2);
    double ram3_(double &ecm, double *xm, double *p1, double *p2, double *p3);
    double ram4_(double &ecm, double *xm, double *p1, double *p2, double *p3, double *p4);
    double ram5_(double &ecm, double *xm, double *p1, double *p2, double *p3, double *p4, double *p5);
}

extern "C" {

    double rndm2_(double dummy) {
        return ((double) rand() / (double) (RAND_MAX));
    }
}

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

void write_tuple_to_file(TNtuple *tup, string var, string file_name, int nBins = 30) {
    double _min = tup->GetMinimum(var.c_str()), _max = tup->GetMaximum(var.c_str());
    TH1F *h = new TH1F("h", "h", nBins, _min, _max);
    tup->Project("h", var.c_str(), "matr2*wt");
    write_histogram_to_file(*h, file_name);
    h->Delete();
}

double const PI = acos(-1), TWO_PI = 2 * PI;

string to_str(int i) {
    char buffer[100];
    sprintf(buffer, "%d", i);
    return buffer;
};

double chi_masses[3] = {Mchi0, Mchi1, Mchi2};
double chi_widths[3] = {10.5e-3, 0.84e-3, 1.93e-3};
double brs_gamma[3] = {1.27e-2, 33.9e-2, 19.2e-2};
double br_mm[3] = {2.2e-4, 5.1e-4, 6.4e-4};

void test_chiGammaJ(int J, int nEv) {
    double pPsi[4], k[4];
    const int nOut = 2;
    double XM[nOut] = {Mpsi, 0};
    double sum = 0;
    double _Mchi = chi_masses[J];
    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = ram2_(_Mchi, XM, pPsi, k);
        wt *= pow(TWO_PI, 4 - 3 * nOut);
        double matr2;
        if (J == 0) matr2 = matr2gamma_0(pPsi, k);
        else if (J == 1) matr2 = matr2gamma_1(pPsi, k);
        else if (J == 2) matr2 = matr2gamma_2(pPsi, k);
        sum += wt*matr2;
    };
    double gamma = sum / ((2 * J + 1)*2 * _Mchi) / nEv;
    double theory = chi_widths[J] * brs_gamma[J];
    //    cout<<"gamma="<<gamma<<endl;
    cout << "chi_c" << J << "-> psi gamma: gamma/theory=" << gamma / theory << endl;

}

void test_chiJ(int J, int nEv, int iDebug = 0) {
    cout << "*********************** chi_c" << J << " ***********************" << endl;
    double theory[3] = {3.95316e-8, 3.77115e-7, 2.705e-7};
    if (J < 0 || J > 2) {
        cout << "Wrong spin J=" << J << endl;
        ::abort();
    };
    string strJ = to_str(J);

    TNtuple tup(("chic" + strJ).c_str(), ("chic" + strJ).c_str(), "q2:m2PsiK1:matr2:wt");
    double p[4], k1[4], k2[4];
    const int nOut = 3;
    double XM[nOut] = {Mpsi, mmu, mmu};
    double sum = 0;
    double _Mchi = chi_masses[J];
    for (int iEv = 0; iEv < nEv; ++iEv) {
        if (iEv % (nEv / 10) == 0) {
            cout << "========= " << (int) (100. * iEv / nEv) << " % ==========" << endl;
        };
        double wt = ram3_(_Mchi, XM, p, k1, k2);
        wt *= pow(TWO_PI, 4 - 3 * nOut);

        double matr2;
        if (J == 0) matr2 = matr2_0(p, k1, k2);
        else if (J == 1) matr2 = matr2_1(p, k1, k2);
        else if (J == 2) matr2 = matr2_2(p, k1, k2);

        double q2 = sum_mass2(k1, k2);
        double m2PsiK1 = sum_mass2(p, k1);
        if (iEv < iDebug) {
            cout << "======== Debug print at iEv=" << iEv << endl;
            println_v4("p", p);
            println_v4("k1", k1);
            println_v4("k2", k2);
            cout << " q2=" << q2 << ";" << endl;
            cout << " m2PsiK1=" << m2PsiK1 << ";" << endl;
            cout << " $$matr2=" << matr2 << ";" << endl;
            cout << "Print[sp[q]/q2];" << endl;
            cout << "Print[sp[p + k1]/m2PsiK1];" << endl;
            cout << "Print[$$matr2/$matr2];" << endl;
        };
        tup.Fill(q2, m2PsiK1, matr2, wt);
        sum += wt*matr2;
    };
    cout << "=============" << endl;
    tup.Write();
    double gamma = sum / ((2 * J + 1)*2 * _Mchi) / nEv;
    double th = theory[J];
    double th_paper = chi_widths[J] * brs_gamma[J] * br_mm[J];
    cout << "chi_c" << J << ": gamma=" << gamma << " vs theoretical " << th << endl;
    cout << "gamma/th=" << gamma / th << endl;
    cout << "gamma/th_paper=" << gamma / th_paper << endl;
    write_tuple_to_file(&tup, "q2", "hQ2_chi" + strJ + "_matr2.hst");
    write_tuple_to_file(&tup, "m2PsiK1", "hm2PsiK1_chi" + strJ + "_matr2.hst");
}

void test_chiJ_psi(int J, int nEv, int iDebug = 0) {
    double theory[3] = {3.95316e-8, 3.77115e-7, 2.705e-7};
    if (J < 0 || J > 1) {
        cout << "Wrong spin J=" << J << endl;
        ::abort();
    };
    string strJ = to_str(J);
    const char* tupName = ("chic" + strJ + "_psi").c_str();

    TNtuple tup(tupName, tupName, "q2:m2PsiK1:m2K1KK1:matr2:wt:wtPsi");

    double p[4], k1[4], k2[4], kk1[4], kk2[4];
    const int nOut = 3, nOutPsi = 2;
    double XM[nOut] = {Mpsi, mmu, mmu}, XMpsi[nOutPsi] = {mmu, mmu};
    double sum = 0;
    double _Mchi = chi_masses[J];
    for (int iEv = 0; iEv < nEv; ++iEv) {
        if (iEv % (nEv / 10) == 0) {
            cout << "========= " << (int) (100. * iEv / nEv) << " % ==========" << endl;
        };
        double wt = ram3_(_Mchi, XM, p, k1, k2);
        wt *= pow(TWO_PI, 4 - 3 * nOut);
        double wtPsi = ram2_(Mpsi, XMpsi, kk1, kk2);
        wtPsi *= pow(TWO_PI, 4 - 3 * nOutPsi);
        apply_boost_to(p, kk1);
        apply_boost_to(p, kk2);

        double matr2;
        if (J == 0) matr2 = matr2_0(kk1, kk2, k1, k2);
        else if (J == 1) matr2 = matr2_1(kk1, kk2, k1, k2);

        double q2 = sum_mass2(k1, k2);
        double m2PsiK1 = sum_mass2(p, k1);
        double m2K1KK1 = sum_mass2(k1, kk1);
        if (iEv < iDebug) {
            cout << "======== Debug print at iEv=" << iEv << endl;
            println_v4("p", p);
            println_v4("k1", k1);
            println_v4("k2", k2);
            println_v4("kk1", kk1);
            println_v4("kk2", kk2);
            cout << " q2=" << q2 << ";" << endl;
            cout << " m2PsiK1=" << m2PsiK1 << ";" << endl;
            cout << " (* m(kk1,kk2)=" << sqrt(sum_mass2(kk1, kk2)) << "*)" << endl;
            cout << " m2K1KK1=" << m2K1KK1 << ";" << endl;
            cout << "Print[sp[q]/q2];" << endl;
            cout << "Print[sp[p + k1]/m2PsiK1];" << endl;
            cout << "Print[sp[k1+kk1]/m2K1KK1];" << endl;
            cout << "Print[" << matr2 << "/($$matr2)];" << endl;
        };
        tup.Fill(q2, m2PsiK1, m2K1KK1, matr2, wt, wtPsi);
        sum += wt*matr2;
    };
    cout << "=============" << endl;
    tup.Write();
    double gamma = sum / ((2 * J + 1)*2 * _Mchi) / nEv;
    double th = theory[J];
    double th_paper = chi_widths[J] * brs_gamma[J] * br_mm[J];
    //    cout << "chi_c" << J << ": gamma=" << gamma << " vs theoretical " << th << endl;
    //    cout << "gamma/th=" << gamma / th << endl;
    //    cout << "gamma/th_paper=" << gamma / th_paper << endl;
    write_tuple_to_file(&tup, "q2", "hQ2_chi" + strJ + "_matr2_psi.hst");
    write_tuple_to_file(&tup, "m2PsiK1", "hm2PsiK1_chi" + strJ + "_matr2_psi.hst");
    write_tuple_to_file(&tup, "m2K1KK1", "hm2K1KK1_chi" + strJ + "_matr2_psi.hst");

}

int main(void) {
    //    test_2body(1e4);
    TFile file("matr2_chic.root", "RECREATE");
    int nEv = 1e6;

    test_chiJ_psi(1, nEv, 10);
    test_chiJ(1, nEv);
    //    test_chiGammaJ(0,nEv);
    //    test_chiGammaJ(1,nEv);
    //    test_chiGammaJ(2,nEv);

    //    test_chiJ(0, nEv);
    //    test_chiJ(1, nEv);
    //    test_chiJ(2, nEv);

    file.Save();
}

