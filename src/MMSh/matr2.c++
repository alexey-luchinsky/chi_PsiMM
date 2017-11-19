#include "kinematics.h"
#include "matr2.h"
#include "math.h"

double mmu = 0.10565837, mmu$2 = pow(mmu, 2), mmu$4 = pow(mmu, 4);
double Mpsi = 3.0969160, Mpsi$2 = pow(Mpsi, 2), Mpsi$4 = pow(Mpsi, 4);
double Mchi0 = 3.4147501, Mchi1 = 3.5106599, Mchi2 = 3.5562000;
double g0 = 0.0743617, g0$2 = pow(g0, 2);
double g1 = 0.990406, g1$2 = pow(g1, 2);
double g2=1.48832, g2$2= pow(g2,2);
double gammaPsi = 92e-6, gammaPsi$2 = pow(gammaPsi, 2);
double alpha = 1. / 137;

// ======== chi_c0 -> psi mu mu =================
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

// ======== chi_c0 -> psi gamma =================
double matr2gamma_0(double pPsi[4], double k[4]) {
    return 2*g0$2*pow(Mchi0,2);
}


// ======== chi_c1 -> psi mu mu =================
double matr2_1(double pPsi[4], double k1[4], double k2[4]) {
    double q2 = sum_mass2(k1, k2);
    double m2PsiK1 = sum_mass2(pPsi, k1);
    double const Mchi = Mchi1, Mchi$2 = pow(Mchi, 2), Mchi$4 = pow(Mchi, 4), Mchi$6 = pow(Mchi, 6);
    double matr2 = 4*alpha*PI*pow(Mchi,-2)*pow(Mpsi,-2)*pow(q2,-2)*(Mchi$6*(2*mmu$2 + q2) - Mchi$4*(2*mmu$2*(Mpsi$2 + q2) + q2*(2*m2PsiK1 - Mpsi$2 + 2*q2)) + 
     Mpsi$2*(2*mmu$4*q2 + 2*mmu$2*(Mpsi$4 - 2*m2PsiK1*q2 - Mpsi$2*q2) + q2*(-2*m2PsiK1*(Mpsi$2 - q2) + 2*pow(m2PsiK1,2) + pow(Mpsi$2 - q2,2))) + 
     Mchi$2*(2*mmu$4*q2 - 2*mmu$2*(Mpsi$4 + 2*m2PsiK1*q2 - 6*Mpsi$2*q2) + q2*(-4*m2PsiK1*Mpsi$2 + Mpsi$4 + 2*m2PsiK1*q2 + 4*Mpsi$2*q2 + 2*pow(m2PsiK1,2) + pow(q2,2))));
    matr2 *= g1$2;
    return matr2;
}

double matr2gamma_1(double pPsi[4], double k[4]) {
    double Mchi=Mchi1, Mchi$2=pow(Mchi,2);
    return g1$2*((Mchi$2 + Mpsi$2)*pow(Mchi,-2)*pow(Mpsi,-2)*pow(Mchi$2 - Mpsi$2,2))/2.;
}

// ======== chi_c2 -> psi mu mu =================
double matr2_2(double pPsi[4], double k1[4], double k2[4]) {
    double q2 = sum_mass2(k1, k2);
    double m2PsiK1 = sum_mass2(pPsi, k1);
    double const Mchi = Mchi2, Mchi$2 = pow(Mchi, 2), Mchi$4 = pow(Mchi, 4), Mchi$6 = pow(Mchi, 6);
    double matr2=(alpha*PI*pow(Mchi,-6)*pow(Mpsi,-2)*pow(q2,-2)*(2*(mmu$2*(28*Mpsi$2 - 9*q2) + (-3*m2PsiK1 + 17*Mpsi$2 - 6*q2)*q2)*pow(Mchi,8) + 3*(2*mmu$2 + q2)*pow(Mchi,10) - 
       8*(m2PsiK1 - mmu$2)*(m2PsiK1 - mmu$2 - Mpsi$2 + q2)*pow(Mpsi$2 - q2,2)*pow(q2,2) + 
       2*Mchi$6*(3*mmu$4*q2 - mmu$2*(62*Mpsi$4 + 6*m2PsiK1*q2 + 19*Mpsi$2*q2 - 5*pow(q2,2)) + 
          q2*(-37*m2PsiK1*Mpsi$2 + 3*Mpsi$4 + 13*m2PsiK1*q2 - 38*Mpsi$2*q2 + 3*pow(m2PsiK1,2) + 9*pow(q2,2))) + 
       2*Mchi$4*(2*mmu$4*(17*Mpsi$2 - 5*q2)*q2 + mmu$2*(105*Mpsi$4*q2 - 2*Mpsi$2*q2*(34*m2PsiK1 + 29*q2) + 28*pow(Mpsi,6) + (20*m2PsiK1 + 9*q2)*pow(q2,2)) + 
          q2*(-8*Mpsi$4*q2 + 2*(17*Mpsi$2 - 5*q2)*pow(m2PsiK1,2) + 17*pow(Mpsi,6) + m2PsiK1*(-37*Mpsi$4 + 86*Mpsi$2*q2 - 21*pow(q2,2)) + 25*Mpsi$2*pow(q2,2) - 
             6*pow(q2,3))) + Mchi$2*(2*mmu$4*q2*(3*Mpsi$4 - 42*Mpsi$2*q2 + 11*pow(q2,2)) + 
          2*mmu$2*(-(Mpsi$4*q2*(6*m2PsiK1 + 55*q2)) + 3*q2*pow(Mpsi,6) + 3*pow(Mpsi,8) + Mpsi$2*(84*m2PsiK1 + 61*q2)*pow(q2,2) - 2*(11*m2PsiK1 + 6*q2)*pow(q2,3)) + 
          q2*(pow(Mpsi$2 - q2,2)*(3*Mpsi$4 - 2*Mpsi$2*q2 + 3*pow(q2,2)) + pow(m2PsiK1,2)*(6*Mpsi$4 - 84*Mpsi$2*q2 + 22*pow(q2,2)) + 
             m2PsiK1*(98*Mpsi$4*q2 - 6*pow(Mpsi,6) - 122*Mpsi$2*pow(q2,2) + 30*pow(q2,3))))))/6.;
    matr2 *= g2$2;
    return matr2;
};

double matr2gamma_2(double pPsi[4], double k[4]) {
    double Mchi=Mchi2, Mchi$2=pow(Mchi,2), Mchi$4=pow(Mchi,4);
    return g2$2*((3*Mchi$4 + 34*Mchi$2*Mpsi$2 + 3*Mpsi$4)*pow(Mchi,-4)*pow(Mpsi,-2)*pow(Mchi$2 - Mpsi$2,2))/48.;
}
