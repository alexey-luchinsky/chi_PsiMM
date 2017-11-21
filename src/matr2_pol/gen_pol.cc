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

using namespace std;

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
        double amp2 = 0;
        for (int iPsi = 0; iPsi < psi->getSpinStates(); ++iPsi) {
            EvtVector4C epsPsi=psi->epsParent(iPsi);
            for (int iMplus = 0; iMplus < mu1->getSpinStates(); ++iMplus) {
                EvtDiracSpinor spMplus = mu1->spParent(iMplus);
                for (int iMminus = 0; iMminus < mu2->getSpinStates(); ++iMminus) {
                    EvtDiracSpinor spMminus = mu2->spParent(iMminus);
                    EvtVector4C epsGamma = EvtLeptonVCurrent(spMminus, spMplus);
                    EvtComplex amp=epsGamma*epsPsi;
                    amp2 += abs2(amp);
                };
            };
        };
        if (iEv < 10) {
            cout << " pPsi=" << pPsi << " m=" << pPsi.mass() << endl;
            cout << " k1=" << k1 << " m=" << k1.mass() << endl;
            cout << " k2=" << k2 << " m=" << k2.mass() << endl;
            Ptot = pPsi + k1 + k2;
            cout << "Ptot=" << Ptot << endl;
            cout << " amp2=" << amp2 << endl;
            double th = 4*(2*(pPsi*k1)*(pPsi*k2)/mPsi/mPsi+k1*k2+3*mmu*mmu);
            cout << " amp2/th=" << amp2 / th << endl;
        };
    }
}
