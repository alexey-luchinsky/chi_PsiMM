#include <iostream>

#include "TFile.h"
#include "TNtuple.h"
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

extern "C" {double rndm2_(double dummy) {
        return ((double) rand() / (double) (RAND_MAX));
    }
}


double const PI = acos(-1), TWO_PI = 2 * PI;

string to_str(int i) {
    char buffer[100];
    sprintf(buffer, "%d", i);
    return buffer;
};

void test_chiJ(int J, int nEv, int iDebug = 0) {
    cout << "*********************** chi_c"<<J<<" ***********************" << endl;
    double chi_masses[3] = {Mchi0, Mchi1, Mchi2};
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
    cout << "chi_c" << J << ": gamma=" << gamma << " vs theoretical " << th << endl;
    cout << "gamma/th=" << gamma / th << endl;


}

int main(void) {
    //    test_2body(1e4);
    TFile file("matr2_chic.root", "RECREATE");
    int nEv = 1e6;
    test_chiJ(0, nEv);
    test_chiJ(1, nEv);
    test_chiJ(2, nEv);
    file.Save();
}

