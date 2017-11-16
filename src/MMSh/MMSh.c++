#include <iostream>
#include "math.h"

#include "TFile.h"
#include "TNtuple.h"

using namespace std;

extern "C" {
    void init_(void);
    double ram2_(double &ecm, double *xm, double *p1, double *p2);
    double ram3_(double &ecm, double *xm, double *p1, double *p2, double *p3);
    double ram4_(double &ecm, double *xm, double *p1, double *p2, double *p3, double *p4);
    double ram5_(double &ecm, double *xm, double *p1, double *p2, double *p3, double *p4, double *p5);
    double rndm2_(double);
}

double rndm2_(double dummy) {
    return ((double) rand() / (double) (RAND_MAX));
}

void print_v4(double p[4]) {
    cout << "{" << p[3] << "," << p[0] << "," << p[1] << "," << p[2] << "}";
}

void print_v4(string name, double p[4]) {
    cout << name << "=";
    print_v4(p);
}

void println_v4(double p[4]) {
    print_v4(p);
    cout << endl;
};

void println_v4(string name, double p[4]) {
    cout << name << "=";
    println_v4(p);
}

double sp(double p1[4], double p2[4]) {
    return p1[3] * p2[3] - p1[0] * p2[0] - p1[1] * p2[1] - p1[2] * p2[2];
}

double sp(double p[4]) {
    return sp(p, p);
}

void sum(double *p1, double *p2, double *P) {
    for (int i = 0; i < 4; ++i) P[i] = p1[i] + p2[i];
}

void sum(double *p1, double *p2, double *p3, double *P) {
    for (int i = 0; i < 4; ++i) P[i] = p1[i] + p2[i] + p3[i];
}

double sum_mass2(double *p1, double *p2) {
    double P[4];
    sum(p1, p2, P);
    return sp(P);
}

double sum_mass2(double *p1, double *p2, double *p3) {
    double P[4];
    sum(p1, p2, p3, P);
    return sp(P);
}
double const PI = acos(-1), TWO_PI = 2 * PI;

void test_2body(int nEv) {
    const int nOut = 2;
    double mmu = 10;
    double XM[nOut] = {0, 0};
    double p1[4], p2[4];

    double sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = ram2_(mmu, XM, p1, p2);
        wt *= pow(TWO_PI, 4 - 3 * nOut);
        sum += wt;
    };
    sum /= nEv;
    cout << " sum=" << sum << " vs " << 1. / (8 * PI) << endl;

}

void test_3body(int nEv) {
    const int nOut = 3;
    double mmu = 10;
    double XM[nOut] = {0, 0, 0};
    double p1[4], p2[4], p3[4], P[4];

    double sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = ram3_(mmu, XM, p1, p2, p3);
        wt *= pow(TWO_PI, 4 - 3 * nOut);
        sum += wt;
    };
    sum /= nEv;
    cout << " sum=" << sum << endl;
}

double mmu = 0.10565837, Mpsi = 3.0969160;
double Mchi0 = 3.4147501, Mchi1 = 3.5106599, Mchi2 = 3.5562000;

double matr2_0(double pPsi[4], double k1[4], double k2[4]) {
    double q2 = sum_mass2(k1, k2);
    double m2PsiK1 = sum_mass2(pPsi, k1);
    double alpha = 1. / 137;
    double const Mchi = Mchi0;
    return 16 * alpha * PI * pow(q2, -2)*(pow(Mchi, 4)*(q2 + 2 * pow(mmu, 2)) + 2 * q2 * pow(mmu, 4) -
            2 * pow(Mchi, 2)*(q2 * (m2PsiK1 + q2) + pow(mmu, 2)*(q2 + 2 * pow(Mpsi, 2))) +
            2 * pow(mmu, 2)*(-2 * m2PsiK1 * q2 + 3 * q2 * pow(Mpsi, 2) + pow(Mpsi, 4)) + q2 * (2 * pow(m2PsiK1, 2) -
            2 * m2PsiK1 * (-q2 + pow(Mpsi, 2)) + pow(q2 + pow(Mpsi, 2), 2))) *
            pow(q2 - pow(Mchi, 2) + pow(Mpsi, 2), -2);
}

double matr2_1(double pPsi[4], double k1[4], double k2[4]) {
    double q2 = sum_mass2(k1, k2);
    double m2PsiK1 = sum_mass2(pPsi, k1);
    double alpha = 1. / 137;
    double const Mchi = Mchi1;
    return 4 * alpha * PI * pow(Mchi, -2) * pow(Mpsi, -2) * pow(q2, -2)*(pow(Mchi, 6)*(q2 + 2 * pow(mmu, 2)) -
            pow(Mchi, 4)*(q2 * (2 * m2PsiK1 + 2 * q2 - pow(Mpsi, 2)) + 2 * pow(mmu, 2)*(q2 + pow(Mpsi, 2))) +
            pow(Mchi, 2)*(2 * q2 * pow(mmu, 4) - 2 * pow(mmu, 2)*(2 * m2PsiK1 * q2 - 6 * q2 * pow(Mpsi, 2) + pow(Mpsi, 4)) +
            q2 * (2 * m2PsiK1 * q2 + 2 * pow(m2PsiK1, 2) - 4 * m2PsiK1 * pow(Mpsi, 2) + 4 * q2 * pow(Mpsi, 2) + pow(Mpsi, 4) + pow(q2, 2))) +
            pow(Mpsi, 2)*(2 * q2 * pow(mmu, 4) + 2 * pow(mmu, 2)*(-2 * m2PsiK1 * q2 - q2 * pow(Mpsi, 2) + pow(Mpsi, 4)) +
            q2 * (2 * pow(m2PsiK1, 2) - 2 * m2PsiK1 * (-q2 + pow(Mpsi, 2)) + pow(-q2 + pow(Mpsi, 2), 2))));
}

double matr2_2(double pPsi[4], double k1[4], double k2[4]) {
    double q2 = sum_mass2(k1, k2);
    double m2PsiK1 = sum_mass2(pPsi, k1);
    double alpha = 1. / 137;
    double const Mchi = Mchi2;
    return (alpha * PI * pow(Mchi, -4) * pow(Mpsi, -2) * pow(q2, -2)*(3 * pow(Mchi, 10)*(q2 + 2 * pow(mmu, 2)) +
            2 * pow(Mchi, 8)*(q2 * (-3 * m2PsiK1 - 6 * q2 + 17 * pow(Mpsi, 2)) + pow(mmu, 2)*(-9 * q2 + 28 * pow(Mpsi, 2))) +
            2 * pow(Mchi, 6)*(3 * q2 * pow(mmu, 4) - pow(mmu, 2)*(6 * m2PsiK1 * q2 + 19 * q2 * pow(Mpsi, 2) + 62 * pow(Mpsi, 4) - 5 * pow(q2, 2)) +
            q2 * (13 * m2PsiK1 * q2 + 3 * pow(m2PsiK1, 2) - 37 * m2PsiK1 * pow(Mpsi, 2) - 38 * q2 * pow(Mpsi, 2) + 3 * pow(Mpsi, 4) + 9 * pow(q2, 2))) +
            2 * pow(Mchi, 4)*(2 * q2 * pow(mmu, 4)*(-5 * q2 + 17 * pow(Mpsi, 2)) + pow(mmu, 2)*
            (-2 * q2 * (34 * m2PsiK1 + 29 * q2) * pow(Mpsi, 2) + 105 * q2 * pow(Mpsi, 4) + 28 * pow(Mpsi, 6) + (20 * m2PsiK1 + 9 * q2) * pow(q2, 2)) +
            q2 * (2 * pow(m2PsiK1, 2)*(-5 * q2 + 17 * pow(Mpsi, 2)) - 8 * q2 * pow(Mpsi, 4) + 17 * pow(Mpsi, 6) + m2PsiK1 * (86 * q2 * pow(Mpsi, 2) - 37 * pow(Mpsi, 4) - 21 * pow(q2, 2)) +
            25 * pow(Mpsi, 2) * pow(q2, 2) - 6 * pow(q2, 3))) - 8 * (m2PsiK1 - pow(mmu, 2))*(m2PsiK1 + q2 - pow(mmu, 2) - pow(Mpsi, 2)) * pow(q2, 2) * pow(-q2 + pow(Mpsi, 2), 2) +
            pow(Mchi, 2)*(2 * q2 * pow(mmu, 4)*(-42 * q2 * pow(Mpsi, 2) + 3 * pow(Mpsi, 4) + 11 * pow(q2, 2)) +
            2 * pow(mmu, 2)*(-(q2 * (6 * m2PsiK1 + 55 * q2) * pow(Mpsi, 4)) + 3 * q2 * pow(Mpsi, 6) + 3 * pow(Mpsi, 8) + (84 * m2PsiK1 + 61 * q2) * pow(Mpsi, 2) * pow(q2, 2) -
            2 * (11 * m2PsiK1 + 6 * q2) * pow(q2, 3)) + q2 * (pow(m2PsiK1, 2)*(-84 * q2 * pow(Mpsi, 2) + 6 * pow(Mpsi, 4) + 22 * pow(q2, 2)) +
            m2PsiK1 * (98 * q2 * pow(Mpsi, 4) - 6 * pow(Mpsi, 6) - 122 * pow(Mpsi, 2) * pow(q2, 2) + 30 * pow(q2, 3)) +
            (-2 * q2 * pow(Mpsi, 2) + 3 * pow(Mpsi, 4) + 3 * pow(q2, 2)) * pow(-q2 + pow(Mpsi, 2), 2))))) / 3.;
}

void test_chi0(int nEv) {
    cout<<"*********************** chi_c0 ***********************"<<endl;
    TNtuple tup("chic0", "chic0", "q2:m2PsiK1:matr2:wt");
    double p[4], k1[4], k2[4];
    const int nOut = 3;
    double XM[nOut] = {Mpsi, mmu, mmu};
    double sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        if (iEv % (nEv / 10) == 0) {
            cout << "========= " << (int) (100. * iEv / nEv) << " % ==========" << endl;
        };
        double wt = ram3_(Mchi0, XM, p, k1, k2);
        double matr2 = matr2_0(p, k1, k2);
        double q2 = sum_mass2(k1, k2);
        double m2PsiK1 = sum_mass2(p, k1);
        tup.Fill(q2,m2PsiK1,matr2,wt);
        sum += wt*matr2;
    };
    sum /= nEv;
    tup.Write();
    cout << "chi_c0: sum=" << sum << endl;
}

void test_chi1(int nEv) {
    cout<<"*********************** chi_c1 ***********************"<<endl;
    TNtuple tup("chic1", "chic1", "q2:m2PsiK1:matr2:wt");
    double p[4], k1[4], k2[4];
    const int nOut = 3;
    double XM[nOut] = {Mpsi, mmu, mmu};
    double sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        if (iEv % (nEv / 10) == 0) {
            cout << "========= " << (int) (100. * iEv / nEv) << " % ==========" << endl;
        };
        double wt = ram3_(Mchi1, XM, p, k1, k2);
        double matr2 = matr2_1(p, k1, k2);
        double q2 = sum_mass2(k1, k2);
        double m2PsiK1 = sum_mass2(p, k1);
        tup.Fill(q2,m2PsiK1,matr2,wt);
        sum += wt*matr2;
    };
    sum /= nEv;
    tup.Write();
    cout << "chi_c0: sum=" << sum << endl;
}

void test_chi2(int nEv) {
    cout<<"*********************** chi_c2 ***********************"<<endl;
    TNtuple tup("chic2", "chic2", "q2:m2PsiK1:matr2:wt");
    double p[4], k1[4], k2[4];
    const int nOut = 3;
    double XM[nOut] = {Mpsi, mmu, mmu};
    double sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        if (iEv % (nEv / 10) == 0) {
            cout << "========= " << (int) (100. * iEv / nEv) << " % ==========" << endl;
        };
        double wt = ram3_(Mchi2, XM, p, k1, k2);
        double matr2 = matr2_2(p, k1, k2);
        double q2 = sum_mass2(k1, k2);
        double m2PsiK1 = sum_mass2(p, k1);
        tup.Fill(q2,m2PsiK1,matr2,wt);
        sum += wt*matr2;
    };
    sum /= nEv;
    tup.Write();
    cout << "chi_c0: sum=" << sum << endl;
}


int main(void) {
    //    test_2body(1e4);
    TFile file("matr2_chic0.root", "RECREATE");
    int nEv=1e6;
    test_chi0(nEv);
    test_chi1(nEv);
    test_chi2(nEv);
    file.Save();
}
