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
    double rndm2_(double);
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




void test_chi2(int nEv) {
    cout << "*********************** chi_c2 ***********************" << endl;
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
        tup.Fill(q2, m2PsiK1, matr2, wt);
        sum += wt*matr2;
    };
    tup.Write();
    sum /= nEv;
}

void test_chi0(int nEv) {
    cout << "*********************** chi_c0 ***********************" << endl;
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
        wt *= pow(TWO_PI, 4 - 3 * nOut);
        double matr2 = matr2_0(p, k1, k2);
        double q2 = sum_mass2(k1, k2);
        double m2PsiK1 = sum_mass2(p, k1);
        if (iEv < 0) {
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
    double gamma = sum / (2 * Mchi0) / nEv;
    double th = 10.5e-3 * 1.26e-2 * 2.97e-4;
    cout << "chi_c0: gamma=" << gamma << " vs theoretical " << th << endl;
    cout << "gamma/th=" << gamma / th << endl;
}


void test_chi1(int nEv) {
    cout << "*********************** chi_c1 ***********************" << endl;
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
        wt *= pow(TWO_PI, 4 - 3 * nOut);
        double matr2 = matr2_1(p, k1, k2);
        double q2 = sum_mass2(k1, k2);
        double m2PsiK1 = sum_mass2(p, k1);
        if (iEv < 10) {
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
    double gamma = sum / ((2*1+1)*2 * Mchi1) / nEv;
    double th = 0.84e-3 * 33.9e-2 * 5.1e-4*2.592;
    cout << "chi_c1: gamma=" << gamma << " vs theoretical " << th << endl;
    cout << "gamma/th=" << gamma / th << endl;

}


int main(void) {
    //    test_2body(1e4);
    TFile file("matr2_chic0.root", "RECREATE");
    int nEv = 1e6;
    test_chi1(nEv);
//    test_chi0_LL(nEv);
    //    test_chi1(nEv);
    //    test_chi2(nEv);
    file.Save();
}

void test_chi0_LL(int nEv) {
    cout << "*********************** chi_c0_LL ***********************" << endl;
    TNtuple tup("chic0LL", "chic0LL", "q2:m2PsiK1:m2K1KK1:matr2:wt");
    double p[4], k1[4], k2[4], kk1[4], kk2[4];
    const int nOut = 3, nOutPsi=2;
    double XM[nOut] = {Mpsi, mmu, mmu};
    double XMpsi[nOutPsi]={mmu,mmu};
    double sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        if (iEv % (nEv / 10) == 0) {
            cout << "========= " << (int) (100. * iEv / nEv) << " % ==========" << endl;
        };
        double wt = ram3_(Mchi0, XM, p, k1, k2);
        wt *= pow(TWO_PI, 4 - 3 * nOut);
        double wtPsi=ram2_(Mpsi,XMpsi,kk1,kk2);
        wtPsi *= pow(TWO_PI, 4 - 3 * nOutPsi);
        apply_boost_to(p,kk1);
        apply_boost_to(p,kk2);
        double matr2=_matr2LL_0(kk1,kk2,k1,k2);
        double q2 = sum_mass2(k1, k2);
        double m2PsiK1 = sum_mass2(p, k1);
        double m2K1KK1=sum_mass2(k1,kk1);
        tup.Fill(q2, m2PsiK1, m2K1KK1, matr2, wt);
        sum += wt*matr2;

        if(iEv<0) {
            cout<<"Debug print at iEv="<<iEv<<endl;
            print_v4("p",p); cout<<"(*m="<<sqrt(sp(p))<<"*)"<<endl;
            print_v4("kk1",kk1); cout<<"(* m="<<sqrt(sp(kk1))<<"*)"<<endl;
            print_v4("kk2",kk2); cout<<" (*m="<<sqrt(sp(kk2))<<"*)"<<endl;
        }
    };
    cout << "=============" << endl;
    tup.Write();
    double gamma = sum / (2 * Mchi0) / nEv;
    double th = 5.971e-2*10.5e-3 * 1.26e-2 * 2.97e-4;
    cout << "chi_c0LL: gamma=" << gamma << " vs theoretical " << th << endl;
    cout << "gamma/th=" << gamma / th << endl;
}
