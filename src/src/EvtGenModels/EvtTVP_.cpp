//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtSVP.cc
//
// Description: Routine to implement radiative decay chi_c0 -> psi gamma
//
//
// Modification history:
//	AVL	Jul 6, 2012	mode created
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtTensorParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"


#include "EvtGenModels/EvtTVP_.hh"


#include <string>
#include <iostream>

using namespace std;

EvtTVP_::~EvtTVP_() {

}

std::string EvtTVP_::getName() {
    return "TVP_";
}

EvtDecayBase* EvtTVP_::clone() {
    return new EvtTVP_;

}

void EvtTVP_::decay_2body(EvtParticle* root) {
    double amp2 = 0;
    root ->initializePhaseSpace(getNDaug(), getDaugs());

    EvtVector4R p = root->getDaug(1)->getP4(), // J/psi momentum
            k = root->getDaug(0)->getP4(); // Photon momentum
    for (int iPsi = 0; iPsi < 3; iPsi++) {
        for (int iGamma = 0; iGamma < 2; iGamma++) {
            for (int iChi = 0; iChi < 5; iChi++) {
                EvtTensor4C epsChi = root->epsTensor(iChi);
                EvtVector4C epsPsi = root->getDaug(1)->epsParent(iPsi).conj();
                EvtVector4C epsGamma = root->getDaug(0)->epsParentPhoton(iGamma).conj();

                // [Baranov, (11)
                // matr = p^mu epsPsi^a epsChi_{a b} ( k_mu epsGamma_b  - k_b epsGamma_mu


                EvtVector4C eee = epsChi.cont1(epsPsi);
                EvtVector4C vvv = (p * k) * eee - (k * eee) * p;
                EvtComplex amp = vvv*epsGamma;
                vertex(iChi, iGamma, iPsi, amp);
                amp2 = amp2 + abs2(amp);
            };
        };
    };
}

void EvtTVP_::decay_3body(EvtParticle* root) {
    root ->initializePhaseSpace(getNDaug(), getDaugs());

    EvtVector4R p = root->getDaug(0)->getP4(), // J/psi momentum
            k1 = root->getDaug(1)->getP4(), // mu+ momentum
            k2 = root->getDaug(2)->getP4(), // mu- momentum
            k = k1 + k2; // photon momentum

    double kSq = k*k;
    if (kSq < 1e-10) {
        return;
    }
    double dSq = delta*delta;
    double dSqDenom = dSq - k.mass2();
    if (fabs(dSqDenom) < 1e-10) {
        return;
    }

    double factor = dSq / (dSqDenom * kSq);

    int iPols[4];
    for (int iChi = 0; iChi < 5; iChi++) {
        iPols[0] = iChi;
        EvtTensor4C epsChi = root->epsTensor(iChi);
        for (int iPsi = 0; iPsi < 3; iPsi++) {
            iPols[1] = iPsi;
            EvtVector4C epsPsi = root->getDaug(0)->epsParent(iPsi).conj();
            for (int iMplus = 0; iMplus < 2; ++iMplus) {
                iPols[2] = iMplus;
                EvtDiracSpinor spMplus = root->getDaug(1)->spParent(iMplus);
                for (int iMminus = 0; iMminus < 2; ++iMminus) {
                    iPols[3] = iMminus;
                    EvtDiracSpinor spMminus = root->getDaug(2)->spParent(iMminus);
                    EvtVector4C epsGamma = EvtLeptonVCurrent(spMplus, spMminus);

                    // [Baranov, (11)
                    // matr = p^mu epsPsi^a epsChi_{a b} ( k_mu epsGamma_b  - k_b epsGamma_mu

                    EvtVector4C eee = epsChi.cont1(epsPsi);
                    EvtVector4C vvv = (p * k) * eee - (k * eee) * p;
                    EvtComplex amp = vvv*epsGamma;
                    amp *= factor;
                    vertex(iPols, amp);
                };
            };
        };
    };
}

void EvtTVP_::decay(EvtParticle *root) {
    if (getNDaug() == 2) decay_2body(root);
    else if (getNDaug() == 3) decay_3body(root);
}

void EvtTVP_::init() {
    checkSpinParent(EvtSpinType::TENSOR);
    if (getNDaug() == 2) { // old SVP
        checkNArg(0);
        checkSpinDaughter(0, EvtSpinType::PHOTON);
        checkSpinDaughter(1, EvtSpinType::VECTOR);
    } else if (getNDaug() == 3) { // chi -> psi mu mu
        checkNDaug(3);
        checkSpinDaughter(0, EvtSpinType::VECTOR);
        checkSpinDaughter(1, EvtSpinType::DIRAC);
        checkSpinDaughter(2, EvtSpinType::DIRAC);
        checkNArg(1);
        delta = getArg(0);
    }

}

void EvtTVP_::initProbMax() {
    if (getNDaug() == 2) setProbMax(2);
    else if (getNDaug() == 3) {
        double ffCor = pow(pow(delta, 2) / (pow(delta, 2) - 0.2), 2);
        if (getDaug(1).getId() == EvtPDL::getId("mu+").getId() && getParentId() == EvtPDL::getId("chi_c2")) {
            setProbMax(ffCor * 520); // tested on 1e6 events
        }
        if (getDaug(1).getId() == EvtPDL::getId("mu+").getId() && getParentId() == EvtPDL::getId("chi_b2")) {
            setProbMax(ffCor * 3000); // tested on 1e6 events
        }
        if (getDaug(1).getId() == EvtPDL::getId("e+").getId()) {
            setProbMax(ffCor * 2e6); // tested on 1e6 events
        };
    };
}


