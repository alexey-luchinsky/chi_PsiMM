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

using namespace std;

int main(void) {
    cout << "gen_pol" << endl;
    int nOut = 3;
    double masses[30];
    masses[0] = 3.1;
    masses[1] = 0.105;
    masses[2] = 0.105;
    double mp = 3.55;
    EvtVector4R p[30];
    EvtRandomEngine *eng = new EvtSimpleRandomEngine();
    EvtRandom::setRandomEngine(eng);
    for (int iEv = 0; iEv < 10; ++iEv) {
        cout << "======== event " << iEv << "=========" << endl;
        EvtGenKine::PhaseSpace(nOut, masses, p, mp);
        EvtVector4R pPsi = p[0], k1 = p[1], k2 = p[2];
        cout << " pPsi=" << pPsi << " m=" << pPsi.mass() << endl;
        cout << " k1=" << k1 << " m=" << k1.mass() << endl;
        cout << " k2=" << k2 << " m=" << k2.mass() << endl;
        EvtVector4R Ptot = pPsi + k1 + k2;
        cout << "Ptot=" << Ptot << endl;
    }
}
