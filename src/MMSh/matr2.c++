#include "kinematics.h"
#include "matr2.h"
#include "math.h"

double mmu = 0.10565837, mmu$2 = pow(mmu, 2), mmu$4 = pow(mmu, 4);
double Mpsi = 3.0969160, Mpsi$2 = pow(Mpsi, 2), Mpsi$4 = pow(Mpsi, 4);
double Mchi0 = 3.4147501, Mchi1 = 3.5106599, Mchi2 = 3.5562000;
double g0 = 0.0743617, g0$2 = pow(g0, 2);
double g1 = 0.990406, g1$2 = pow(g1, 2);
double g2 = 1.48832, g2$2 = pow(g2, 2);
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
    return 2 * g0$2 * pow(Mchi0, 2);
}

// ======== chi_c0 -> psi mu mu -> mu(kk1) mu(kk2) mu(k1) mu(k2)
double matr2_0_mm(double kk1[4], double kk2[4], double k1[4], double k2[4]) {
    double Mchi=Mchi0, Mchi$2=pow(Mchi,2);
    double p[4]; sum(kk1,kk2,p);
    double m2PsiK1=sum_mass2(p,k1), m2PsiK2=sum_mass2(p,k2), q2=sum_mass2(k1,k2);
    double k2kk1=sp(k2,kk1),  //
            k2kk2=sp(k2,kk2), //
            kk1kk2=sp(kk1,kk2), //
            kk1p=sp(kk1,p), kk2p=sp(kk2,p), k1kk2=sp(k1,kk2),
            k1kk1=sp(k1,kk1);
//    cout<<"(*======================*)"<<endl;
//    println_v4("kk1",kk1);
//    println_v4("kk2",kk2);
//    println_v4("p",p);
//    println_v4("k1",k1);
//    println_v4("k2",k2);
//    cout<<"Print["<<m2PsiK1<<"/(sp[p+k1])];"<<endl;
//    cout<<"Print["<<m2PsiK2<<"/(sp[p+k2])];"<<endl;
//    cout<<"Print["<<q2<<"/(sp[k1+k2])];"<<endl;
//    cout<<"Print["<<k2kk1<<"/(sp[k2,kk1])];"<<endl;
//    cout<<"Print["<<k2kk2<<"/(sp[k2,kk2])];"<<endl;
//    cout<<"Print["<<kk1kk2<<"/(sp[kk1,kk2])];"<<endl;
//    cout<<"Print["<<kk1p<<"/(sp[kk1,p])];"<<endl;
//    cout<<"Print["<<kk2p<<"/(sp[kk2,p])];"<<endl;
//    cout<<"Print["<<k1kk2<<"/(sp[k1,kk2])];"<<endl;
//    cout<<"Print["<<k1kk1<<"/(sp[k1,kk1])];"<<endl;


    double matr2=
64*alpha*Mchi$2*PI*pow(m2PsiK1 + m2PsiK2 - 2*(mmu$2 + Mpsi$2),-2)*pow(q2,-2)*
   (8*k2kk1*k2kk2*m2PsiK1*mmu$2 + 4*kk1kk2*m2PsiK1*m2PsiK2*mmu$2 - 4*k2kk1*k2kk2*mmu$4 - 8*kk1kk2*m2PsiK1*mmu$4 - 8*kk1kk2*m2PsiK2*mmu$4 + 4*m2PsiK1*m2PsiK2*mmu$4 + 
     8*k2kk1*k2kk2*m2PsiK1*Mpsi$2 - 8*k2kk1*k2kk2*mmu$2*Mpsi$2 - 8*kk1kk2*m2PsiK1*mmu$2*Mpsi$2 - 8*kk1kk2*m2PsiK2*mmu$2*Mpsi$2 + 16*kk1kk2*mmu$4*Mpsi$2 - 
     8*m2PsiK1*mmu$4*Mpsi$2 - 8*m2PsiK2*mmu$4*Mpsi$2 - 4*k2kk1*k2kk2*Mpsi$4 + 8*kk1kk2*mmu$2*Mpsi$4 + 8*mmu$4*Mpsi$4 + 2*k2kk2*kk1p*m2PsiK1*q2 + 
     2*k2kk1*kk2p*m2PsiK1*q2 + 2*k2kk2*kk1p*m2PsiK2*q2 + 2*k2kk1*kk2p*m2PsiK2*q2 - 2*kk1kk2*m2PsiK1*m2PsiK2*q2 - 4*k2kk2*kk1p*mmu$2*q2 - 4*k2kk1*kk2p*mmu$2*q2 + 
     2*kk1kk2*m2PsiK1*mmu$2*q2 + 2*kk1kk2*m2PsiK2*mmu$2*q2 - 2*kk1kk2*mmu$4*q2 - 2*m2PsiK1*mmu$4*q2 - 2*m2PsiK2*mmu$4*q2 - 4*k2kk1*k2kk2*Mpsi$2*q2 - 
     4*k2kk2*kk1p*Mpsi$2*q2 - 4*k2kk1*kk2p*Mpsi$2*q2 + 2*kk1kk2*m2PsiK1*Mpsi$2*q2 + 2*kk1kk2*m2PsiK2*Mpsi$2*q2 - 4*kk1kk2*mmu$2*Mpsi$2*q2 - 
     2*m2PsiK1*mmu$2*Mpsi$2*q2 - 2*m2PsiK2*mmu$2*Mpsi$2*q2 + 4*mmu$4*Mpsi$2*q2 - 2*kk1kk2*Mpsi$4*q2 + 2*mmu$2*Mpsi$4*q2 + 
     2*k1kk2*(kk1p*(m2PsiK1 + m2PsiK2 - 2*(mmu$2 + Mpsi$2))*q2 + 2*k2kk1*
         (mmu$4 + m2PsiK1*(m2PsiK2 - mmu$2 - Mpsi$2) + 2*mmu$2*Mpsi$2 - m2PsiK2*(mmu$2 + Mpsi$2) + Mpsi$4 - Mpsi$2*q2)) - 4*k2kk1*k2kk2*pow(m2PsiK1,2) + 
     2*kk1kk2*mmu$2*pow(m2PsiK1,2) + 2*mmu$4*pow(m2PsiK1,2) + mmu$2*q2*pow(m2PsiK1,2) + 2*kk1kk2*mmu$2*pow(m2PsiK2,2) + 2*mmu$4*pow(m2PsiK2,2) + 
     mmu$2*q2*pow(m2PsiK2,2) - 2*k1kk1*(-(kk2p*(m2PsiK1 + m2PsiK2 - 2*(mmu$2 + Mpsi$2))*q2) + 
        2*k2kk2*(-mmu$4 - 2*mmu$2*Mpsi$2 + m2PsiK2*(mmu$2 + Mpsi$2) + m2PsiK1*(-m2PsiK2 + mmu$2 + Mpsi$2) - Mpsi$4 + Mpsi$2*q2) + 
        2*k1kk2*(mmu$4 + 2*mmu$2*Mpsi$2 - 2*m2PsiK2*(mmu$2 + Mpsi$2) + Mpsi$4 + Mpsi$2*q2 + pow(m2PsiK2,2))) + 8*kk1kk2*pow(mmu,6) - 8*m2PsiK1*pow(mmu,6) - 
     8*m2PsiK2*pow(mmu,6) + 16*Mpsi$2*pow(mmu,6) + 2*q2*pow(mmu,6) + 8*pow(mmu,8) + 2*kk1kk2*Mpsi$2*pow(q2,2) + 2*mmu$2*Mpsi$2*pow(q2,2));
    matr2 *= g0$2;
    matr2 /= (pow(sqrt(sp(p))-Mpsi,2)+gammaPsi$2);
//    cout<<"Print["<<matr2<<"/($$matr2)];"<<endl;
    return matr2;
}

// chi_c0 -> mu+ mu- mu+ mu- with symmetrization
double matr2_0(double kk1[4], double kk2[4], double k1[4], double k2[4]) {
    return matr2_0_mm(kk1,kk2,k1,k2)+matr2_0_mm(k1,kk2,kk1,k2)+
            matr2_0_mm(kk1,k2,k1,kk2)+matr2_0_mm(k1,k2,kk1,kk2);
}

// ======== chi_c1 -> psi mu mu =================

double matr2_1(double pPsi[4], double k1[4], double k2[4]) {
    double q2 = sum_mass2(k1, k2);
    double m2PsiK1 = sum_mass2(pPsi, k1);
    double const Mchi = Mchi1, Mchi$2 = pow(Mchi, 2), Mchi$4 = pow(Mchi, 4), Mchi$6 = pow(Mchi, 6);
    double matr2 = 4 * alpha * PI * pow(Mchi, -2) * pow(Mpsi, -2) * pow(q2, -2)*(Mchi$6 * (2 * mmu$2 + q2) - Mchi$4 * (2 * mmu$2 * (Mpsi$2 + q2) + q2 * (2 * m2PsiK1 - Mpsi$2 + 2 * q2)) +
            Mpsi$2 * (2 * mmu$4 * q2 + 2 * mmu$2 * (Mpsi$4 - 2 * m2PsiK1 * q2 - Mpsi$2 * q2) + q2 * (-2 * m2PsiK1 * (Mpsi$2 - q2) + 2 * pow(m2PsiK1, 2) + pow(Mpsi$2 - q2, 2))) +
            Mchi$2 * (2 * mmu$4 * q2 - 2 * mmu$2 * (Mpsi$4 + 2 * m2PsiK1 * q2 - 6 * Mpsi$2 * q2) + q2 * (-4 * m2PsiK1 * Mpsi$2 + Mpsi$4 + 2 * m2PsiK1 * q2 + 4 * Mpsi$2 * q2 + 2 * pow(m2PsiK1, 2) + pow(q2, 2))));
    matr2 *= g1$2;
    return matr2;
}

double matr2gamma_1(double pPsi[4], double k[4]) {
    double Mchi = Mchi1, Mchi$2 = pow(Mchi, 2);
    return g1$2 * ((Mchi$2 + Mpsi$2) * pow(Mchi, -2) * pow(Mpsi, -2) * pow(Mchi$2 - Mpsi$2, 2)) / 2.;
}

// ======== chi_c1 -> psi mu mu -> mu(kk1) mu(kk2) mu(k1) mu(k2)
double matr2_1_mm(double kk1[4], double kk2[4], double k1[4], double k2[4]) {
    double Mchi=Mchi1, Mchi$2=pow(Mchi,2), Mchi$4=pow(Mchi,4);
    double p[4], P[4]; sum(kk1,kk2,p); sum(p,k1,k2,P);
    double m2PsiK1=sum_mass2(p,k1), m2PsiK2=sum_mass2(p,k2), q2=sum_mass2(k1,k2);
    double k2kk1=sp(k2,kk1),  //
            k2kk2=sp(k2,kk2), //
            kk1kk2=sp(kk1,kk2), //
            kk1p=sp(kk1,p), kk2p=sp(kk2,p), k1kk2=sp(k1,kk2),
            k1kk1=sp(k1,kk1),
            kk1P=sp(kk1,P), kk2P=sp(kk2,P);
    

    double matr2=
16*alpha*PI*pow(Mchi,-2)*pow(q2,-2)*(-8*k2kk1*k2kk2*m2PsiK2*Mchi$2 + 4*k2kk1*k2kk2*Mchi$4 - 8*k2kk2*kk1P*m2PsiK1*mmu$2 - 8*k2kk1*kk2P*m2PsiK1*mmu$2 - 
     8*k2kk1*k2kk2*m2PsiK2*mmu$2 - 8*k2kk2*kk1P*m2PsiK2*mmu$2 - 8*k2kk1*kk2P*m2PsiK2*mmu$2 - 4*kk1kk2*m2PsiK1*m2PsiK2*mmu$2 + 8*k2kk1*k2kk2*Mchi$2*mmu$2 + 
     16*k2kk2*kk1P*Mchi$2*mmu$2 + 16*k2kk1*kk2P*Mchi$2*mmu$2 + 8*kk1kk2*m2PsiK1*Mchi$2*mmu$2 + 8*kk1kk2*m2PsiK2*Mchi$2*mmu$2 - 8*kk1kk2*Mchi$4*mmu$2 + 
     4*k2kk1*k2kk2*mmu$4 + 16*k2kk2*kk1P*mmu$4 + 16*k2kk1*kk2P*mmu$4 + 8*kk1kk2*m2PsiK1*mmu$4 + 8*kk1kk2*m2PsiK2*mmu$4 + 4*m2PsiK1*m2PsiK2*mmu$4 - 
     16*kk1kk2*Mchi$2*mmu$4 - 8*m2PsiK1*Mchi$2*mmu$4 - 8*m2PsiK2*Mchi$2*mmu$4 + 8*Mchi$4*mmu$4 - 2*k2kk2*kk1P*m2PsiK1*q2 - 2*k2kk1*kk2P*m2PsiK1*q2 + 
     2*k2kk2*kk1P*m2PsiK2*q2 + 2*k2kk1*kk2P*m2PsiK2*q2 + 2*kk1kk2*m2PsiK1*m2PsiK2*q2 + 4*k2kk1*k2kk2*Mchi$2*q2 - 2*kk1kk2*m2PsiK1*Mchi$2*q2 - 
     2*kk1kk2*m2PsiK2*Mchi$2*q2 + 2*kk1kk2*Mchi$4*q2 - 16*kk1P*kk2P*mmu$2*q2 - 2*kk1kk2*m2PsiK1*mmu$2*q2 - 2*kk1kk2*m2PsiK2*mmu$2*q2 + 12*kk1kk2*Mchi$2*mmu$2*q2 - 
     2*m2PsiK1*Mchi$2*mmu$2*q2 - 2*m2PsiK2*Mchi$2*mmu$2*q2 + 2*Mchi$4*mmu$2*q2 + 2*kk1kk2*mmu$4*q2 - 2*m2PsiK1*mmu$4*q2 - 2*m2PsiK2*mmu$4*q2 + 12*Mchi$2*mmu$4*q2 - 
     2*k1kk2*(2*k2kk1*(Mchi$4 + m2PsiK1*(m2PsiK2 - Mchi$2 - mmu$2) + 2*Mchi$2*mmu$2 - m2PsiK2*(Mchi$2 + mmu$2) + mmu$4 - Mchi$2*q2) + 
        kk1P*(-8*mmu$2*(Mchi$2 + mmu$2) + m2PsiK1*(4*mmu$2 - q2) + m2PsiK2*(4*mmu$2 + q2))) - 2*kk1kk2*mmu$2*pow(m2PsiK1,2) + 2*mmu$4*pow(m2PsiK1,2) + 
     mmu$2*q2*pow(m2PsiK1,2) + 2*k1kk1*(2*k2kk2*(-Mchi$4 - 2*Mchi$2*mmu$2 + m2PsiK2*(Mchi$2 + mmu$2) + m2PsiK1*(-m2PsiK2 + Mchi$2 + mmu$2) - mmu$4 + Mchi$2*q2) + 
        kk2P*(8*mmu$2*(Mchi$2 + mmu$2) + m2PsiK1*(-4*mmu$2 + q2) - m2PsiK2*(4*mmu$2 + q2)) + 
        2*k1kk2*(Mchi$4 - 2*m2PsiK1*(Mchi$2 + mmu$2) + mmu$4 + Mchi$2*(2*mmu$2 + q2) + pow(m2PsiK1,2))) + 4*k2kk1*k2kk2*pow(m2PsiK2,2) - 
     2*kk1kk2*mmu$2*pow(m2PsiK2,2) + 2*mmu$4*pow(m2PsiK2,2) + mmu$2*q2*pow(m2PsiK2,2) - 8*kk1kk2*pow(mmu,6) - 8*m2PsiK1*pow(mmu,6) - 8*m2PsiK2*pow(mmu,6) + 
     16*Mchi$2*pow(mmu,6) + 2*q2*pow(mmu,6) + 8*pow(mmu,8) + 2*kk1kk2*Mchi$2*pow(q2,2) + 6*Mchi$2*mmu$2*pow(q2,2));
    matr2 *= g1$2;
    matr2 /= (pow(sqrt(sp(p))-Mpsi,2)+gammaPsi$2);
//    cout<<"Print["<<matr2<<"/($$matr2)];"<<endl;
    return matr2;
}

// chi_c0 -> mu+ mu- mu+ mu- with symmetrization
double matr2_1(double kk1[4], double kk2[4], double k1[4], double k2[4]) {
    return matr2_1_mm(kk1,kk2,k1,k2)+matr2_1_mm(k1,kk2,kk1,k2)+
            matr2_1_mm(kk1,k2,k1,kk2)+matr2_1_mm(k1,k2,kk1,kk2);
}


// ======== chi_c2 -> psi mu mu =================

double matr2_2(double pPsi[4], double k1[4], double k2[4]) {
    double q2 = sum_mass2(k1, k2);
    double m2PsiK1 = sum_mass2(pPsi, k1);
    double const Mchi = Mchi2, Mchi$2 = pow(Mchi, 2), Mchi$4 = pow(Mchi, 4), Mchi$6 = pow(Mchi, 6);
    double matr2 = (alpha * PI * pow(Mchi, -6) * pow(Mpsi, -2) * pow(q2, -2)*(2 * (mmu$2 * (28 * Mpsi$2 - 9 * q2) + (-3 * m2PsiK1 + 17 * Mpsi$2 - 6 * q2) * q2) * pow(Mchi, 8) + 3 * (2 * mmu$2 + q2) * pow(Mchi, 10) -
            8 * (m2PsiK1 - mmu$2)*(m2PsiK1 - mmu$2 - Mpsi$2 + q2) * pow(Mpsi$2 - q2, 2) * pow(q2, 2) +
            2 * Mchi$6 * (3 * mmu$4 * q2 - mmu$2 * (62 * Mpsi$4 + 6 * m2PsiK1 * q2 + 19 * Mpsi$2 * q2 - 5 * pow(q2, 2)) +
            q2 * (-37 * m2PsiK1 * Mpsi$2 + 3 * Mpsi$4 + 13 * m2PsiK1 * q2 - 38 * Mpsi$2 * q2 + 3 * pow(m2PsiK1, 2) + 9 * pow(q2, 2))) +
            2 * Mchi$4 * (2 * mmu$4 * (17 * Mpsi$2 - 5 * q2) * q2 + mmu$2 * (105 * Mpsi$4 * q2 - 2 * Mpsi$2 * q2 * (34 * m2PsiK1 + 29 * q2) + 28 * pow(Mpsi, 6) + (20 * m2PsiK1 + 9 * q2) * pow(q2, 2)) +
            q2 * (-8 * Mpsi$4 * q2 + 2 * (17 * Mpsi$2 - 5 * q2) * pow(m2PsiK1, 2) + 17 * pow(Mpsi, 6) + m2PsiK1 * (-37 * Mpsi$4 + 86 * Mpsi$2 * q2 - 21 * pow(q2, 2)) + 25 * Mpsi$2 * pow(q2, 2) -
            6 * pow(q2, 3))) + Mchi$2 * (2 * mmu$4 * q2 * (3 * Mpsi$4 - 42 * Mpsi$2 * q2 + 11 * pow(q2, 2)) +
            2 * mmu$2 * (-(Mpsi$4 * q2 * (6 * m2PsiK1 + 55 * q2)) + 3 * q2 * pow(Mpsi, 6) + 3 * pow(Mpsi, 8) + Mpsi$2 * (84 * m2PsiK1 + 61 * q2) * pow(q2, 2) - 2 * (11 * m2PsiK1 + 6 * q2) * pow(q2, 3)) +
            q2 * (pow(Mpsi$2 - q2, 2)*(3 * Mpsi$4 - 2 * Mpsi$2 * q2 + 3 * pow(q2, 2)) + pow(m2PsiK1, 2)*(6 * Mpsi$4 - 84 * Mpsi$2 * q2 + 22 * pow(q2, 2)) +
            m2PsiK1 * (98 * Mpsi$4 * q2 - 6 * pow(Mpsi, 6) - 122 * Mpsi$2 * pow(q2, 2) + 30 * pow(q2, 3)))))) / 6.;
    matr2 *= g2$2;
    return matr2;
};

double matr2gamma_2(double pPsi[4], double k[4]) {
    double Mchi = Mchi2, Mchi$2 = pow(Mchi, 2), Mchi$4 = pow(Mchi, 4);
    return g2$2 * ((3 * Mchi$4 + 34 * Mchi$2 * Mpsi$2 + 3 * Mpsi$4) * pow(Mchi, -4) * pow(Mpsi, -2) * pow(Mchi$2 - Mpsi$2, 2)) / 48.;
}

// ======== chi_c2 -> psi mu mu -> mu(kk1) mu(kk2) mu(k1) mu(k2)
double matr2_2_mm(double kk1[4], double kk2[4], double k1[4], double k2[4]) {
    double Mchi=Mchi2, Mchi$2=pow(Mchi,2), Mchi$4=pow(Mchi,4), Mchi$6=pow(Mchi,6);
    double p[4], P[4]; sum(kk1,kk2,p); sum(p,k1,k2,P);
    double m2PsiK1=sum_mass2(p,k1), m2PsiK2=sum_mass2(p,k2), q2=sum_mass2(k1,k2);
    double k2kk1=sp(k2,kk1),  //
            k2kk2=sp(k2,kk2), //
            kk1kk2=sp(kk1,kk2), //
            kk1p=sp(kk1,p), kk2p=sp(kk2,p), k1kk2=sp(k1,kk2),
            k1kk1=sp(k1,kk1),
            kk1P=sp(kk1,P), kk2P=sp(kk2,P);
    

    double matr2=
(-8*alpha*PI*pow(Mchi,-6)*pow(q2,-2)*(2*k2kk2*kk1P*m2PsiK1*m2PsiK2*Mchi$4 + 2*k2kk1*kk2P*m2PsiK1*m2PsiK2*Mchi$4 - 8*kk1P*kk2P*m2PsiK1*m2PsiK2*Mchi$4 - 
       4*kk1kk2*m2PsiK1*m2PsiK2*Mchi$6 + 4*k2kk2*kk1P*m2PsiK1*m2PsiK2*Mchi$2*mmu$2 + 4*k2kk1*kk2P*m2PsiK1*m2PsiK2*Mchi$2*mmu$2 - 
       56*kk1P*kk2P*m2PsiK1*m2PsiK2*Mchi$2*mmu$2 - 8*k2kk1*k2kk2*m2PsiK1*Mchi$4*mmu$2 + 2*k2kk2*kk1P*m2PsiK1*Mchi$4*mmu$2 + 2*k2kk1*kk2P*m2PsiK1*Mchi$4*mmu$2 - 
       2*k2kk2*kk1P*m2PsiK2*Mchi$4*mmu$2 - 2*k2kk1*kk2P*m2PsiK2*Mchi$4*mmu$2 - 32*kk1kk2*m2PsiK1*m2PsiK2*Mchi$4*mmu$2 - 10*m2PsiK1*m2PsiK2*Mchi$6*mmu$2 - 
       32*kk1P*kk2P*m2PsiK1*m2PsiK2*mmu$4 + 4*k2kk2*kk1P*m2PsiK1*Mchi$2*mmu$4 + 4*k2kk1*kk2P*m2PsiK1*Mchi$2*mmu$4 + 48*kk1P*kk2P*m2PsiK1*Mchi$2*mmu$4 - 
       4*k2kk2*kk1P*m2PsiK2*Mchi$2*mmu$4 - 4*k2kk1*kk2P*m2PsiK2*Mchi$2*mmu$4 + 48*kk1P*kk2P*m2PsiK2*Mchi$2*mmu$4 - 16*kk1kk2*m2PsiK1*m2PsiK2*Mchi$2*mmu$4 + 
       4*k2kk1*k2kk2*Mchi$4*mmu$4 + 32*kk1kk2*m2PsiK1*Mchi$4*mmu$4 + 32*kk1kk2*m2PsiK2*Mchi$4*mmu$4 - 80*m2PsiK1*m2PsiK2*Mchi$4*mmu$4 + 
       2*k2kk2*kk1P*m2PsiK1*m2PsiK2*Mchi$2*Mpsi$2 + 2*k2kk1*kk2P*m2PsiK1*m2PsiK2*Mchi$2*Mpsi$2 - 16*kk1P*kk2P*m2PsiK1*m2PsiK2*Mchi$2*Mpsi$2 - 
       8*k2kk1*k2kk2*m2PsiK1*Mchi$4*Mpsi$2 + 2*k2kk2*kk1P*m2PsiK1*Mchi$4*Mpsi$2 + 2*k2kk1*kk2P*m2PsiK1*Mchi$4*Mpsi$2 - 2*k2kk2*kk1P*m2PsiK2*Mchi$4*Mpsi$2 - 
       2*k2kk1*kk2P*m2PsiK2*Mchi$4*Mpsi$2 - 8*kk1kk2*m2PsiK1*m2PsiK2*Mchi$4*Mpsi$2 - 32*kk1P*kk2P*m2PsiK1*m2PsiK2*mmu$2*Mpsi$2 + 
       6*k2kk2*kk1P*m2PsiK1*Mchi$2*mmu$2*Mpsi$2 + 6*k2kk1*kk2P*m2PsiK1*Mchi$2*mmu$2*Mpsi$2 + 48*kk1P*kk2P*m2PsiK1*Mchi$2*mmu$2*Mpsi$2 - 
       6*k2kk2*kk1P*m2PsiK2*Mchi$2*mmu$2*Mpsi$2 - 6*k2kk1*kk2P*m2PsiK2*Mchi$2*mmu$2*Mpsi$2 + 48*kk1P*kk2P*m2PsiK2*Mchi$2*mmu$2*Mpsi$2 - 
       16*kk1kk2*m2PsiK1*m2PsiK2*Mchi$2*mmu$2*Mpsi$2 + 8*k2kk1*k2kk2*Mchi$4*mmu$2*Mpsi$2 + 32*kk1kk2*m2PsiK1*Mchi$4*mmu$2*Mpsi$2 + 
       32*kk1kk2*m2PsiK2*Mchi$4*mmu$2*Mpsi$2 - 20*m2PsiK1*m2PsiK2*Mchi$4*mmu$2*Mpsi$2 - 96*kk1P*kk2P*Mchi$2*mmu$4*Mpsi$2 - 40*m2PsiK1*m2PsiK2*Mchi$2*mmu$4*Mpsi$2 - 
       64*kk1kk2*Mchi$4*mmu$4*Mpsi$2 + 80*m2PsiK1*Mchi$4*mmu$4*Mpsi$2 + 80*m2PsiK2*Mchi$4*mmu$4*Mpsi$2 - 8*kk1P*kk2P*m2PsiK1*m2PsiK2*Mpsi$4 + 
       2*k2kk2*kk1P*m2PsiK1*Mchi$2*Mpsi$4 + 2*k2kk1*kk2P*m2PsiK1*Mchi$2*Mpsi$4 - 2*k2kk2*kk1P*m2PsiK2*Mchi$2*Mpsi$4 - 2*k2kk1*kk2P*m2PsiK2*Mchi$2*Mpsi$4 - 
       4*kk1kk2*m2PsiK1*m2PsiK2*Mchi$2*Mpsi$4 + 4*k2kk1*k2kk2*Mchi$4*Mpsi$4 - 48*kk1P*kk2P*Mchi$2*mmu$2*Mpsi$4 - 10*m2PsiK1*m2PsiK2*Mchi$2*mmu$2*Mpsi$4 - 
       32*kk1kk2*Mchi$4*mmu$2*Mpsi$4 - 80*Mchi$4*mmu$4*Mpsi$4 - 2*kk1P*kk2p*m2PsiK1*m2PsiK2*Mchi$2*q2 - 2*kk1p*kk2P*m2PsiK1*m2PsiK2*Mchi$2*q2 + 
       12*kk1P*kk2P*m2PsiK1*m2PsiK2*Mchi$2*q2 - 2*k2kk2*kk1p*m2PsiK1*Mchi$4*q2 + k2kk2*kk1P*m2PsiK1*Mchi$4*q2 - 2*k2kk1*kk2p*m2PsiK1*Mchi$4*q2 + 
       2*kk1P*kk2p*m2PsiK1*Mchi$4*q2 + k2kk1*kk2P*m2PsiK1*Mchi$4*q2 + 2*kk1p*kk2P*m2PsiK1*Mchi$4*q2 - 8*kk1P*kk2P*m2PsiK1*Mchi$4*q2 - 2*k2kk2*kk1p*m2PsiK2*Mchi$4*q2 + 
       k2kk2*kk1P*m2PsiK2*Mchi$4*q2 - 2*k2kk1*kk2p*m2PsiK2*Mchi$4*q2 + 2*kk1P*kk2p*m2PsiK2*Mchi$4*q2 + k2kk1*kk2P*m2PsiK2*Mchi$4*q2 + 2*kk1p*kk2P*m2PsiK2*Mchi$4*q2 - 
       8*kk1P*kk2P*m2PsiK2*Mchi$4*q2 + 10*kk1kk2*m2PsiK1*m2PsiK2*Mchi$4*q2 - 4*kk1kk2*m2PsiK1*Mchi$6*q2 - 4*kk1kk2*m2PsiK2*Mchi$6*q2 + 
       4*kk1P*kk2p*m2PsiK1*Mchi$2*mmu$2*q2 + 4*kk1p*kk2P*m2PsiK1*Mchi$2*mmu$2*q2 - 12*kk1P*kk2P*m2PsiK1*Mchi$2*mmu$2*q2 + 4*kk1P*kk2p*m2PsiK2*Mchi$2*mmu$2*q2 + 
       4*kk1p*kk2P*m2PsiK2*Mchi$2*mmu$2*q2 - 12*kk1P*kk2P*m2PsiK2*Mchi$2*mmu$2*q2 + 4*k2kk2*kk1p*Mchi$4*mmu$2*q2 - 2*k2kk2*kk1P*Mchi$4*mmu$2*q2 + 
       4*k2kk1*kk2p*Mchi$4*mmu$2*q2 - 4*kk1P*kk2p*Mchi$4*mmu$2*q2 - 2*k2kk1*kk2P*Mchi$4*mmu$2*q2 - 4*kk1p*kk2P*Mchi$4*mmu$2*q2 + 16*kk1P*kk2P*Mchi$4*mmu$2*q2 - 
       12*kk1kk2*m2PsiK1*Mchi$4*mmu$2*q2 - 12*kk1kk2*m2PsiK2*Mchi$4*mmu$2*q2 + 20*m2PsiK1*m2PsiK2*Mchi$4*mmu$2*q2 + 8*kk1kk2*Mchi$6*mmu$2*q2 - 
       10*m2PsiK1*Mchi$6*mmu$2*q2 - 10*m2PsiK2*Mchi$6*mmu$2*q2 - 4*kk1P*kk2p*Mchi$2*mmu$4*q2 - 4*kk1p*kk2P*Mchi$2*mmu$4*q2 + 12*kk1P*kk2P*Mchi$2*mmu$4*q2 + 
       12*kk1kk2*Mchi$4*mmu$4*q2 - 20*m2PsiK1*Mchi$4*mmu$4*q2 - 20*m2PsiK2*Mchi$4*mmu$4*q2 + 20*Mchi$6*mmu$4*q2 + 16*kk1P*kk2P*m2PsiK1*m2PsiK2*Mpsi$2*q2 + 
       3*k2kk2*kk1P*m2PsiK1*Mchi$2*Mpsi$2*q2 + 2*kk1P*kk2p*m2PsiK1*Mchi$2*Mpsi$2*q2 + 3*k2kk1*kk2P*m2PsiK1*Mchi$2*Mpsi$2*q2 + 2*kk1p*kk2P*m2PsiK1*Mchi$2*Mpsi$2*q2 - 
       28*kk1P*kk2P*m2PsiK1*Mchi$2*Mpsi$2*q2 + 3*k2kk2*kk1P*m2PsiK2*Mchi$2*Mpsi$2*q2 + 2*kk1P*kk2p*m2PsiK2*Mchi$2*Mpsi$2*q2 + 3*k2kk1*kk2P*m2PsiK2*Mchi$2*Mpsi$2*q2 + 
       2*kk1p*kk2P*m2PsiK2*Mchi$2*Mpsi$2*q2 - 28*kk1P*kk2P*m2PsiK2*Mchi$2*Mpsi$2*q2 + 8*kk1kk2*m2PsiK1*m2PsiK2*Mchi$2*Mpsi$2*q2 + 4*k2kk1*k2kk2*Mchi$4*Mpsi$2*q2 + 
       4*k2kk2*kk1p*Mchi$4*Mpsi$2*q2 - 6*k2kk2*kk1P*Mchi$4*Mpsi$2*q2 + 4*k2kk1*kk2p*Mchi$4*Mpsi$2*q2 - 4*kk1P*kk2p*Mchi$4*Mpsi$2*q2 - 6*k2kk1*kk2P*Mchi$4*Mpsi$2*q2 - 
       4*kk1p*kk2P*Mchi$4*Mpsi$2*q2 + 32*kk1P*kk2P*Mchi$4*Mpsi$2*q2 - 20*kk1kk2*m2PsiK1*Mchi$4*Mpsi$2*q2 - 20*kk1kk2*m2PsiK2*Mchi$4*Mpsi$2*q2 + 
       16*kk1kk2*Mchi$6*Mpsi$2*q2 - 32*kk1P*kk2P*m2PsiK1*mmu$2*Mpsi$2*q2 - 32*kk1P*kk2P*m2PsiK2*mmu$2*Mpsi$2*q2 - 6*k2kk2*kk1P*Mchi$2*mmu$2*Mpsi$2*q2 - 
       4*kk1P*kk2p*Mchi$2*mmu$2*Mpsi$2*q2 - 6*k2kk1*kk2P*Mchi$2*mmu$2*Mpsi$2*q2 - 4*kk1p*kk2P*Mchi$2*mmu$2*Mpsi$2*q2 + 56*kk1P*kk2P*Mchi$2*mmu$2*Mpsi$2*q2 - 
       16*kk1kk2*m2PsiK1*Mchi$2*mmu$2*Mpsi$2*q2 - 16*kk1kk2*m2PsiK2*Mchi$2*mmu$2*Mpsi$2*q2 + 20*m2PsiK1*m2PsiK2*Mchi$2*mmu$2*Mpsi$2*q2 + 
       40*kk1kk2*Mchi$4*mmu$2*Mpsi$2*q2 - 40*m2PsiK1*Mchi$4*mmu$2*Mpsi$2*q2 - 40*m2PsiK2*Mchi$4*mmu$2*Mpsi$2*q2 + 40*Mchi$6*mmu$2*Mpsi$2*q2 + 
       32*kk1P*kk2P*mmu$4*Mpsi$2*q2 + 16*kk1kk2*Mchi$2*mmu$4*Mpsi$2*q2 - 40*m2PsiK1*Mchi$2*mmu$4*Mpsi$2*q2 - 40*m2PsiK2*Mchi$2*mmu$4*Mpsi$2*q2 + 
       80*Mchi$4*mmu$4*Mpsi$2*q2 - 8*kk1P*kk2P*m2PsiK1*Mpsi$4*q2 - 8*kk1P*kk2P*m2PsiK2*Mpsi$4*q2 - 2*k2kk2*kk1P*Mchi$2*Mpsi$4*q2 - 2*k2kk1*kk2P*Mchi$2*Mpsi$4*q2 + 
       12*kk1P*kk2P*Mchi$2*Mpsi$4*q2 - 4*kk1kk2*m2PsiK1*Mchi$2*Mpsi$4*q2 - 4*kk1kk2*m2PsiK2*Mchi$2*Mpsi$4*q2 + 12*kk1kk2*Mchi$4*Mpsi$4*q2 + 
       16*kk1P*kk2P*mmu$2*Mpsi$4*q2 + 8*kk1kk2*Mchi$2*mmu$2*Mpsi$4*q2 - 10*m2PsiK1*Mchi$2*mmu$2*Mpsi$4*q2 - 10*m2PsiK2*Mchi$2*mmu$2*Mpsi$4*q2 + 
       20*Mchi$4*mmu$2*Mpsi$4*q2 + 20*Mchi$2*mmu$4*Mpsi$4*q2 + 8*kk1P*kk2P*m2PsiK2*Mchi$2*pow(m2PsiK1,2) + 4*k2kk1*k2kk2*Mchi$4*pow(m2PsiK1,2) - 
       2*k2kk2*kk1P*Mchi$4*pow(m2PsiK1,2) - 2*k2kk1*kk2P*Mchi$4*pow(m2PsiK1,2) + 4*kk1P*kk2P*Mchi$4*pow(m2PsiK1,2) + 4*kk1kk2*m2PsiK2*Mchi$4*pow(m2PsiK1,2) + 
       2*kk1kk2*Mchi$6*pow(m2PsiK1,2) + 16*kk1P*kk2P*m2PsiK2*mmu$2*pow(m2PsiK1,2) - 6*k2kk2*kk1P*Mchi$2*mmu$2*pow(m2PsiK1,2) - 
       6*k2kk1*kk2P*Mchi$2*mmu$2*pow(m2PsiK1,2) + 4*kk1P*kk2P*Mchi$2*mmu$2*pow(m2PsiK1,2) + 8*kk1kk2*m2PsiK2*Mchi$2*mmu$2*pow(m2PsiK1,2) + 
       10*m2PsiK2*Mchi$4*mmu$2*pow(m2PsiK1,2) + 5*Mchi$6*mmu$2*pow(m2PsiK1,2) + 16*kk1P*kk2P*mmu$4*pow(m2PsiK1,2) + 8*kk1kk2*Mchi$2*mmu$4*pow(m2PsiK1,2) + 
       20*m2PsiK2*Mchi$2*mmu$4*pow(m2PsiK1,2) + 8*kk1P*kk2P*m2PsiK2*Mpsi$2*pow(m2PsiK1,2) - 4*k2kk2*kk1P*Mchi$2*Mpsi$2*pow(m2PsiK1,2) - 
       4*k2kk1*kk2P*Mchi$2*Mpsi$2*pow(m2PsiK1,2) + 8*kk1P*kk2P*Mchi$2*Mpsi$2*pow(m2PsiK1,2) + 4*kk1kk2*m2PsiK2*Mchi$2*Mpsi$2*pow(m2PsiK1,2) + 
       4*kk1kk2*Mchi$4*Mpsi$2*pow(m2PsiK1,2) + 16*kk1P*kk2P*mmu$2*Mpsi$2*pow(m2PsiK1,2) + 8*kk1kk2*Mchi$2*mmu$2*Mpsi$2*pow(m2PsiK1,2) + 
       10*m2PsiK2*Mchi$2*mmu$2*Mpsi$2*pow(m2PsiK1,2) + 10*Mchi$4*mmu$2*Mpsi$2*pow(m2PsiK1,2) + 20*Mchi$2*mmu$4*Mpsi$2*pow(m2PsiK1,2) + 
       4*kk1P*kk2P*Mpsi$4*pow(m2PsiK1,2) + 2*kk1kk2*Mchi$2*Mpsi$4*pow(m2PsiK1,2) + 5*Mchi$2*mmu$2*Mpsi$4*pow(m2PsiK1,2) - kk1P*kk2p*Mchi$2*q2*pow(m2PsiK1,2) - 
       kk1p*kk2P*Mchi$2*q2*pow(m2PsiK1,2) + kk1kk2*Mchi$4*q2*pow(m2PsiK1,2) + 8*kk1P*kk2P*Mpsi$2*q2*pow(m2PsiK1,2) + 4*kk1kk2*Mchi$2*Mpsi$2*q2*pow(m2PsiK1,2) + 
       10*Mchi$2*mmu$2*Mpsi$2*q2*pow(m2PsiK1,2) + 2*k2kk2*kk1P*Mchi$2*pow(m2PsiK1,3) + 2*k2kk1*kk2P*Mchi$2*pow(m2PsiK1,3) - 8*kk1P*kk2P*Mchi$2*pow(m2PsiK1,3) - 
       4*kk1kk2*Mchi$4*pow(m2PsiK1,3) - 16*kk1P*kk2P*mmu$2*pow(m2PsiK1,3) - 8*kk1kk2*Mchi$2*mmu$2*pow(m2PsiK1,3) - 10*Mchi$4*mmu$2*pow(m2PsiK1,3) - 
       20*Mchi$2*mmu$4*pow(m2PsiK1,3) - 8*kk1P*kk2P*Mpsi$2*pow(m2PsiK1,3) - 4*kk1kk2*Mchi$2*Mpsi$2*pow(m2PsiK1,3) - 10*Mchi$2*mmu$2*Mpsi$2*pow(m2PsiK1,3) + 
       4*kk1P*kk2P*pow(m2PsiK1,4) + 2*kk1kk2*Mchi$2*pow(m2PsiK1,4) + 5*Mchi$2*mmu$2*pow(m2PsiK1,4) - 2*k2kk2*kk1P*m2PsiK1*Mchi$2*pow(m2PsiK2,2) - 
       2*k2kk1*kk2P*m2PsiK1*Mchi$2*pow(m2PsiK2,2) + 8*kk1P*kk2P*m2PsiK1*Mchi$2*pow(m2PsiK2,2) + 4*kk1P*kk2P*Mchi$4*pow(m2PsiK2,2) + 
       4*kk1kk2*m2PsiK1*Mchi$4*pow(m2PsiK2,2) + 2*kk1kk2*Mchi$6*pow(m2PsiK2,2) + 16*kk1P*kk2P*m2PsiK1*mmu$2*pow(m2PsiK2,2) + 
       2*k2kk2*kk1P*Mchi$2*mmu$2*pow(m2PsiK2,2) + 2*k2kk1*kk2P*Mchi$2*mmu$2*pow(m2PsiK2,2) + 4*kk1P*kk2P*Mchi$2*mmu$2*pow(m2PsiK2,2) + 
       8*kk1kk2*m2PsiK1*Mchi$2*mmu$2*pow(m2PsiK2,2) + 10*m2PsiK1*Mchi$4*mmu$2*pow(m2PsiK2,2) + 5*Mchi$6*mmu$2*pow(m2PsiK2,2) + 16*kk1P*kk2P*mmu$4*pow(m2PsiK2,2) + 
       8*kk1kk2*Mchi$2*mmu$4*pow(m2PsiK2,2) + 20*m2PsiK1*Mchi$2*mmu$4*pow(m2PsiK2,2) + 8*kk1P*kk2P*m2PsiK1*Mpsi$2*pow(m2PsiK2,2) + 
       2*k2kk2*kk1P*Mchi$2*Mpsi$2*pow(m2PsiK2,2) + 2*k2kk1*kk2P*Mchi$2*Mpsi$2*pow(m2PsiK2,2) + 8*kk1P*kk2P*Mchi$2*Mpsi$2*pow(m2PsiK2,2) + 
       4*kk1kk2*m2PsiK1*Mchi$2*Mpsi$2*pow(m2PsiK2,2) + 4*kk1kk2*Mchi$4*Mpsi$2*pow(m2PsiK2,2) + 16*kk1P*kk2P*mmu$2*Mpsi$2*pow(m2PsiK2,2) + 
       8*kk1kk2*Mchi$2*mmu$2*Mpsi$2*pow(m2PsiK2,2) + 10*m2PsiK1*Mchi$2*mmu$2*Mpsi$2*pow(m2PsiK2,2) + 10*Mchi$4*mmu$2*Mpsi$2*pow(m2PsiK2,2) + 
       20*Mchi$2*mmu$4*Mpsi$2*pow(m2PsiK2,2) + 4*kk1P*kk2P*Mpsi$4*pow(m2PsiK2,2) + 2*kk1kk2*Mchi$2*Mpsi$4*pow(m2PsiK2,2) + 5*Mchi$2*mmu$2*Mpsi$4*pow(m2PsiK2,2) - 
       kk1P*kk2p*Mchi$2*q2*pow(m2PsiK2,2) - kk1p*kk2P*Mchi$2*q2*pow(m2PsiK2,2) + kk1kk2*Mchi$4*q2*pow(m2PsiK2,2) + 8*kk1P*kk2P*Mpsi$2*q2*pow(m2PsiK2,2) + 
       4*kk1kk2*Mchi$2*Mpsi$2*q2*pow(m2PsiK2,2) + 10*Mchi$2*mmu$2*Mpsi$2*q2*pow(m2PsiK2,2) - 8*kk1P*kk2P*pow(m2PsiK1,2)*pow(m2PsiK2,2) - 
       4*kk1kk2*Mchi$2*pow(m2PsiK1,2)*pow(m2PsiK2,2) - 10*Mchi$2*mmu$2*pow(m2PsiK1,2)*pow(m2PsiK2,2) - 8*kk1P*kk2P*Mchi$2*pow(m2PsiK2,3) - 
       4*kk1kk2*Mchi$4*pow(m2PsiK2,3) - 16*kk1P*kk2P*mmu$2*pow(m2PsiK2,3) - 8*kk1kk2*Mchi$2*mmu$2*pow(m2PsiK2,3) - 10*Mchi$4*mmu$2*pow(m2PsiK2,3) - 
       20*Mchi$2*mmu$4*pow(m2PsiK2,3) - 8*kk1P*kk2P*Mpsi$2*pow(m2PsiK2,3) - 4*kk1kk2*Mchi$2*Mpsi$2*pow(m2PsiK2,3) - 10*Mchi$2*mmu$2*Mpsi$2*pow(m2PsiK2,3) + 
       4*kk1P*kk2P*pow(m2PsiK2,4) + 2*kk1kk2*Mchi$2*pow(m2PsiK2,4) + 5*Mchi$2*mmu$2*pow(m2PsiK2,4) - 48*kk1P*kk2P*Mchi$2*pow(mmu,6) - 
       40*m2PsiK1*m2PsiK2*Mchi$2*pow(mmu,6) - 32*kk1kk2*Mchi$4*pow(mmu,6) + 80*m2PsiK1*Mchi$4*pow(mmu,6) + 80*m2PsiK2*Mchi$4*pow(mmu,6) - 
       160*Mchi$4*Mpsi$2*pow(mmu,6) + 20*Mchi$4*q2*pow(mmu,6) + 40*Mchi$2*Mpsi$2*q2*pow(mmu,6) + 20*Mchi$2*pow(m2PsiK1,2)*pow(mmu,6) + 
       20*Mchi$2*pow(m2PsiK2,2)*pow(mmu,6) - 80*Mchi$4*pow(mmu,8) + 
       k1kk2*Mchi$2*(-2*kk1p*Mchi$2*(m2PsiK1 + m2PsiK2 - 2*(mmu$2 + Mpsi$2))*q2 + 
          4*k2kk1*Mchi$2*(-mmu$4 - 2*mmu$2*Mpsi$2 + m2PsiK2*(mmu$2 + Mpsi$2) + m2PsiK1*(-m2PsiK2 + mmu$2 + Mpsi$2) - Mpsi$4 + Mpsi$2*q2) + 
          kk1P*(-2*q2*(Mchi$2*(mmu$2 + 3*Mpsi$2) + Mpsi$4 + mmu$2*(3*Mpsi$2 - q2) - Mpsi$2*q2) + 2*(-m2PsiK2 + mmu$2 + Mpsi$2)*pow(m2PsiK1,2) - 
             2*(Mchi$2 + 3*mmu$2 + 2*Mpsi$2)*pow(m2PsiK2,2) + 2*pow(m2PsiK2,3) + 
             m2PsiK1*(-4*mmu$4 - 6*mmu$2*Mpsi$2 + 2*m2PsiK2*(Mchi$2 + 2*mmu$2 + Mpsi$2) - 2*Mpsi$4 + 3*Mpsi$2*q2 + Mchi$2*(-2*mmu$2 - 2*Mpsi$2 + q2) - pow(q2,2)) + 
             m2PsiK2*(4*mmu$4 + 6*mmu$2*Mpsi$2 + 2*Mpsi$4 + 3*Mpsi$2*q2 + Mchi$2*(2*mmu$2 + 2*Mpsi$2 + q2) - pow(q2,2)))) - 8*kk1P*kk2P*m2PsiK1*m2PsiK2*pow(q2,2) - 
       k2kk2*kk1P*m2PsiK1*Mchi$2*pow(q2,2) - k2kk1*kk2P*m2PsiK1*Mchi$2*pow(q2,2) + 8*kk1P*kk2P*m2PsiK1*Mchi$2*pow(q2,2) - k2kk2*kk1P*m2PsiK2*Mchi$2*pow(q2,2) - 
       k2kk1*kk2P*m2PsiK2*Mchi$2*pow(q2,2) + 8*kk1P*kk2P*m2PsiK2*Mchi$2*pow(q2,2) - 4*kk1kk2*m2PsiK1*m2PsiK2*Mchi$2*pow(q2,2) + 4*kk1kk2*m2PsiK1*Mchi$4*pow(q2,2) + 
       4*kk1kk2*m2PsiK2*Mchi$4*pow(q2,2) + 16*kk1P*kk2P*m2PsiK1*mmu$2*pow(q2,2) + 16*kk1P*kk2P*m2PsiK2*mmu$2*pow(q2,2) + 2*k2kk2*kk1P*Mchi$2*mmu$2*pow(q2,2) + 
       2*k2kk1*kk2P*Mchi$2*mmu$2*pow(q2,2) - 16*kk1P*kk2P*Mchi$2*mmu$2*pow(q2,2) + 8*kk1kk2*m2PsiK1*Mchi$2*mmu$2*pow(q2,2) + 8*kk1kk2*m2PsiK2*Mchi$2*mmu$2*pow(q2,2) - 
       10*m2PsiK1*m2PsiK2*Mchi$2*mmu$2*pow(q2,2) - 8*kk1kk2*Mchi$4*mmu$2*pow(q2,2) + 10*m2PsiK1*Mchi$4*mmu$2*pow(q2,2) + 10*m2PsiK2*Mchi$4*mmu$2*pow(q2,2) - 
       16*kk1P*kk2P*mmu$4*pow(q2,2) - 8*kk1kk2*Mchi$2*mmu$4*pow(q2,2) + 20*m2PsiK1*Mchi$2*mmu$4*pow(q2,2) + 20*m2PsiK2*Mchi$2*mmu$4*pow(q2,2) - 
       20*Mchi$4*mmu$4*pow(q2,2) + 8*kk1P*kk2P*m2PsiK1*Mpsi$2*pow(q2,2) + 8*kk1P*kk2P*m2PsiK2*Mpsi$2*pow(q2,2) + 2*k2kk2*kk1P*Mchi$2*Mpsi$2*pow(q2,2) + 
       2*k2kk1*kk2P*Mchi$2*Mpsi$2*pow(q2,2) - 28*kk1P*kk2P*Mchi$2*Mpsi$2*pow(q2,2) + 4*kk1kk2*m2PsiK1*Mchi$2*Mpsi$2*pow(q2,2) + 
       4*kk1kk2*m2PsiK2*Mchi$2*Mpsi$2*pow(q2,2) - 16*kk1kk2*Mchi$4*Mpsi$2*pow(q2,2) - 16*kk1P*kk2P*mmu$2*Mpsi$2*pow(q2,2) - 8*kk1kk2*Mchi$2*mmu$2*Mpsi$2*pow(q2,2) + 
       10*m2PsiK1*Mchi$2*mmu$2*Mpsi$2*pow(q2,2) + 10*m2PsiK2*Mchi$2*mmu$2*Mpsi$2*pow(q2,2) - 40*Mchi$4*mmu$2*Mpsi$2*pow(q2,2) - 20*Mchi$2*mmu$4*Mpsi$2*pow(q2,2) - 
       4*kk1P*kk2P*pow(m2PsiK1,2)*pow(q2,2) - 2*kk1kk2*Mchi$2*pow(m2PsiK1,2)*pow(q2,2) - 5*Mchi$2*mmu$2*pow(m2PsiK1,2)*pow(q2,2) - 
       4*kk1P*kk2P*pow(m2PsiK2,2)*pow(q2,2) - 2*kk1kk2*Mchi$2*pow(m2PsiK2,2)*pow(q2,2) - 5*Mchi$2*mmu$2*pow(m2PsiK2,2)*pow(q2,2) - 20*Mchi$2*pow(mmu,6)*pow(q2,2) + 
       k1kk1*Mchi$2*(2*kk2P*m2PsiK1*m2PsiK2*Mchi$2 + 4*kk2P*m2PsiK1*m2PsiK2*mmu$2 - 2*kk2P*m2PsiK1*Mchi$2*mmu$2 + 2*kk2P*m2PsiK2*Mchi$2*mmu$2 - 4*kk2P*m2PsiK1*mmu$4 + 
          4*kk2P*m2PsiK2*mmu$4 + 2*kk2P*m2PsiK1*m2PsiK2*Mpsi$2 - 2*kk2P*m2PsiK1*Mchi$2*Mpsi$2 + 2*kk2P*m2PsiK2*Mchi$2*Mpsi$2 - 6*kk2P*m2PsiK1*mmu$2*Mpsi$2 + 
          6*kk2P*m2PsiK2*mmu$2*Mpsi$2 - 2*kk2P*m2PsiK1*Mpsi$4 + 2*kk2P*m2PsiK2*Mpsi$4 - 2*kk2p*m2PsiK1*Mchi$2*q2 + kk2P*m2PsiK1*Mchi$2*q2 - 2*kk2p*m2PsiK2*Mchi$2*q2 + 
          kk2P*m2PsiK2*Mchi$2*q2 + 4*kk2p*Mchi$2*mmu$2*q2 - 2*kk2P*Mchi$2*mmu$2*q2 + 3*kk2P*m2PsiK1*Mpsi$2*q2 + 3*kk2P*m2PsiK2*Mpsi$2*q2 + 4*kk2p*Mchi$2*Mpsi$2*q2 - 
          6*kk2P*Mchi$2*Mpsi$2*q2 - 6*kk2P*mmu$2*Mpsi$2*q2 - 2*kk2P*Mpsi$4*q2 + 
          4*k2kk2*Mchi$2*(-mmu$4 - 2*mmu$2*Mpsi$2 + m2PsiK2*(mmu$2 + Mpsi$2) + m2PsiK1*(-m2PsiK2 + mmu$2 + Mpsi$2) - Mpsi$4 + Mpsi$2*q2) - 
          2*kk2P*m2PsiK2*pow(m2PsiK1,2) + 2*kk2P*mmu$2*pow(m2PsiK1,2) + 2*kk2P*Mpsi$2*pow(m2PsiK1,2) - 2*kk2P*Mchi$2*pow(m2PsiK2,2) - 6*kk2P*mmu$2*pow(m2PsiK2,2) - 
          4*kk2P*Mpsi$2*pow(m2PsiK2,2) + 4*k1kk2*Mchi$2*(mmu$4 + 2*mmu$2*Mpsi$2 - 2*m2PsiK2*(mmu$2 + Mpsi$2) + Mpsi$4 + Mpsi$2*q2 + pow(m2PsiK2,2)) + 
          2*kk2P*pow(m2PsiK2,3) - kk2P*m2PsiK1*pow(q2,2) - kk2P*m2PsiK2*pow(q2,2) + 2*kk2P*mmu$2*pow(q2,2) + 2*kk2P*Mpsi$2*pow(q2,2))))/3.;
    matr2 *= g2$2;
    matr2 /= (pow(sqrt(sp(p))-Mpsi,2)+gammaPsi$2);
//    cout<<"Print["<<matr2<<"/($$matr2)];"<<endl;
    return matr2;
}

// chi_c2 -> mu+ mu- mu+ mu- with symmetrization
double matr2_2(double kk1[4], double kk2[4], double k1[4], double k2[4]) {
    return matr2_2_mm(kk1,kk2,k1,k2)+matr2_2_mm(k1,kk2,kk1,k2)+
            matr2_2_mm(kk1,k2,k1,kk2)+matr2_2_mm(k1,k2,kk1,kk2);
}
