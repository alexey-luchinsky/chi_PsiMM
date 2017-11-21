#include <iostream>
#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtSimpleRandomEngine.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtVectorParticle.hh"
#include "EvtGenBase/EvtDiracParticle.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtId.hh"

#include "../MMSh/matr2.h"
#include "../MMSh/kinematics.h"
using namespace std;

double get_amp2_simple(EvtVectorParticle *psi, EvtDiracParticle *mu1, EvtDiracParticle *mu2) {
    double amp2 = 0;
    for (int iPsi = 0; iPsi < psi->getSpinStates(); ++iPsi) {
        EvtVector4C epsPsi = psi->epsParent(iPsi);
        for (int iMplus = 0; iMplus < mu1->getSpinStates(); ++iMplus) {
            EvtDiracSpinor spMplus = mu1->spParent(iMplus);
            for (int iMminus = 0; iMminus < mu2->getSpinStates(); ++iMminus) {
                EvtDiracSpinor spMminus = mu2->spParent(iMminus);
                EvtVector4C epsGamma = EvtLeptonVCurrent(spMminus, spMplus);
                EvtComplex amp = epsGamma*epsPsi;
                amp2 += abs2(amp);
            };
        };
    };
    return amp2;
}


double get_amp2_0(EvtVectorParticle *psi, EvtDiracParticle *mu1, EvtDiracParticle *mu2) {
    EvtVector4R p = psi->getP4(), // J/psi momentum
            k1 = mu1->getP4(), // mu+ momentum
            k2 = mu2->getP4(), // mu- momentum
            k = k1 + k2; // photon momentum

    double kSq = k*k;
    if (kSq < 1e-10) {
        return 0;
    }
    double kp = k*p;
    if (fabs(kp) < 1e-10) {
        return 0;
    }


    double factor = 1. / kSq;
    double amp2 = 0;
    for (int iPsi = 0; iPsi < psi->getSpinStates(); ++iPsi) {
        EvtVector4C epsPsi = psi->epsParent(iPsi).conj();
        for (int iMplus = 0; iMplus < mu1->getSpinStates(); ++iMplus) {
            EvtDiracSpinor spMplus = mu1->spParent(iMplus);
            for (int iMminus = 0; iMminus < mu2->getSpinStates(); ++iMminus) {
                EvtDiracSpinor spMminus = mu2->spParent(iMminus);
                EvtVector4C epsGamma = EvtLeptonVCurrent(spMplus, spMminus);
                EvtComplex amp = (epsPsi * epsGamma) - (epsPsi * k)*(epsGamma * p) / kp;
                amp *= factor;
                amp2 += abs2(amp);
            };
        };
    };
    return amp2;
}



double matr2_0(EvtVector4R pPsi, EvtVector4R k1, EvtVector4R k2) {
    double q2 = (k1+k2).mass2();
    double m2PsiK1 = (pPsi+k1).mass2();
    double const Mchi = 3.4147501;
    double const Mchi$2 = pow(Mchi, 2);
    double Mpsi = 3.0969160, Mpsi$2 = pow(Mpsi, 2), Mpsi$4 = pow(Mpsi, 4);
    double mmu = 0.10565837, mmu$2 = pow(mmu, 2), mmu$4 = pow(mmu, 4);
    double g0 = 0.0743617, g0$2 = pow(g0, 2);
    double alpha = 1. / 137;
    double PI=acos(-1);

    double matr2 = 16 * alpha * Mchi$2 * PI * pow(q2, -2)*(2 * mmu$4 * q2 + 2 * mmu$2 * (Mpsi$4 - 2 * m2PsiK1 * q2 + 3 * Mpsi$2 * q2) -
            2 * Mchi$2 * (q2 * (m2PsiK1 + q2) + mmu$2 * (2 * Mpsi$2 + q2)) +
            (2 * mmu$2 + q2) * pow(Mchi, 4) + q2 * (-2 * m2PsiK1 * (Mpsi$2 - q2) + 2 * pow(m2PsiK1, 2) +
            pow(Mpsi$2 + q2, 2))) * pow(-Mchi$2 + Mpsi$2 + q2, -2);
    matr2 *= g0$2;
    return matr2;

}

int main(void) {
    cout << "gen_pol" << endl;
    EvtRandomEngine *eng = new EvtSimpleRandomEngine();
    EvtRandom::setRandomEngine(eng);
    EvtGen *myGenerator = new EvtGen("../src/chic_matr2.dec", "../src/evt.pdl", eng);
    int nOut = 3;
    double masses[30];

    double mPsi = EvtPDL::getMass(EvtPDL::getId("J/psi")), mmu = EvtPDL::getMass(EvtPDL::getId("mu+"));
    masses[0] = mPsi;
    masses[1] = mmu;
    masses[2] = mmu;
    double mp = 3.55;
    EvtVector4R p[30];


    EvtScalarParticle *chi0 = new EvtScalarParticle();
    EvtVectorParticle *psi = new EvtVectorParticle();
    EvtDiracParticle *mu1 = new EvtDiracParticle(), *mu2 = new EvtDiracParticle();
    EvtVector4R pPsi, k1, k2, Ptot;
    for (int iEv = 0; iEv < 10; ++iEv) {
        cout << "======== event " << iEv << "=========" << endl;
        EvtGenKine::PhaseSpace(nOut, masses, p, mp);
        pPsi = p[0];
        psi->init(EvtPDL::getId("J/psi"), pPsi);
        k1 = p[1];
        mu1->init(EvtPDL::getId("mu+"), k1);
        k2 = p[2];
        mu2->init(EvtPDL::getId("mu-"), k2);
//        double amp2=get_amp2_simple(psi,mu1,mu2);
        double amp2=get_amp2_0(psi,mu1,mu2);
        if (iEv < 10) {
            cout << " pPsi=" << pPsi << " m=" << pPsi.mass() << endl;
            cout << " k1=" << k1 << " m=" << k1.mass() << endl;
            cout << " k2=" << k2 << " m=" << k2.mass() << endl;
            Ptot = pPsi + k1 + k2;
            cout << "Ptot=" << Ptot << endl;
            cout << " amp2=" << amp2 << endl;
//            double th = 4 * (2 * (pPsi * k1)*(pPsi * k2) / mPsi / mPsi + k1 * k2 + 3 * mmu * mmu);
            double th=matr2_0(pPsi, k1, k2);
            cout << " amp2/th=" << amp2 / th << endl;
        };
    }
}
