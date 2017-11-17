#include "kinematics.h"
#include "matr2.h"
#include "math.h"

double mmu = 0.10565837, mmu$2 = pow(mmu, 2), mmu$4 = pow(mmu, 4);
double Mpsi = 3.0969160, Mpsi$2 = pow(Mpsi, 2), Mpsi$4 = pow(Mpsi, 4);
double Mchi0 = 3.4147501, Mchi1 = 3.5106599, Mchi2 = 3.5562000;
double g0 = 0.0743617, g0$2 = pow(g0, 2);
double gammaPsi=92e-6, gammaPsi$2=pow(gammaPsi,2);
double alpha = 1. / 137;

// ======== chi_c0 =================

double matr2_0(double pPsi[4], double k1[4], double k2[4]) {
    double q2 = sum_mass2(k1, k2);
    double m2PsiK1 = sum_mass2(pPsi, k1);
    double const Mchi = Mchi0;
    double const Mchi$2 = pow(Mchi, 2);
    double matr2 = 16 * alpha * Mchi$2 * PI * pow(q2, -2)*(2 * mmu$4 * q2 + 2 * mmu$2 * (Mpsi$4 - 2 * m2PsiK1 * q2 + 3 * Mpsi$2 * q2) -
            2 * Mchi$2 * (q2 * (m2PsiK1 + q2) + mmu$2 * (2 * Mpsi$2 + q2)) +
            (2 * mmu$2 + q2) * pow(Mchi, 4) + q2 * (-2 * m2PsiK1 * (Mpsi$2 - q2) + 2 * pow(m2PsiK1, 2) +
            pow(Mpsi$2 + q2, 2))) * pow(-Mchi$2 + Mpsi$2 + q2, -2);
    matr2 *= g0$2;
    return matr2;
}

double _matr2LL_0(double kk1[4], double kk2[4], double k1[4], double k2[4]) {
    double q2 = sum_mass2(k1, k2);
    double m2PsiK1 = sum_mass2(kk1,kk2, k1);
    double m2PsiK2 = sum_mass2(kk1,kk2, k2);
    double const Mchi = Mchi0;
    double const Mchi$2 = pow(Mchi, 2);
    double k1kk2=sp(k1,kk2);
    double k2kk1=sp(k2,kk1);
    double k1kk1=sp(k1,kk1);
    double k2kk2=sp(k2,kk2);
//    cout<<"mass(kk1,kk2)="<<sqrt(sum_mass2(kk1,kk2))<<endl;
    double matr2 = 64 * alpha * Mchi$2 * PI * (2 * k1kk2 * k2kk1 + 2 * k1kk1 * k2kk2 +
            mmu$2 * (Mpsi$2 + q2)) * pow(m2PsiK1, 2) * pow(m2PsiK1 + m2PsiK2 -
            2 * (mmu$2 + Mpsi$2), -2) * pow(q2, -2);
    matr2 /= ( pow(sum_mass2(kk1,kk2)-Mpsi$2,2) + gammaPsi$2);
    matr2 /= 1.75636e+14;
    return matr2;
}

double matr2_1(double pPsi[4], double k1[4], double k2[4]) {
    double q2 = sum_mass2(k1, k2);
    double m2PsiK1 = sum_mass2(pPsi, k1);
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

