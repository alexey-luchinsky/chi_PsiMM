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


#include "EvtGenModels/EvtSVP_.hh"


#include <string>
#include <iostream>

using namespace std;

EvtSVP_::~EvtSVP_() {

}

std::string EvtSVP_::getName() {
    return "SVP_";
}

EvtDecayBase* EvtSVP_::clone() {
    return new EvtSVP_;

}

void EvtSVP_::decay_2body(EvtParticle* root) {
    root ->initializePhaseSpace(getNDaug(), getDaugs());

    EvtVector4R p = root->getDaug(1)->getP4(), // J/psi momentum
            k = root->getDaug(0)->getP4(); // Photon momentum
    for (int iPsi = 0; iPsi < 3; iPsi++) {
        for (int iGamma = 0; iGamma < 2; iGamma++) {
            EvtVector4C epsPsi = root->getDaug(1)->epsParent(iPsi).conj();
            EvtVector4C epsGamma = root->getDaug(0)->epsParentPhoton(iGamma).conj();
            EvtComplex amp = (epsPsi * epsGamma) - (epsPsi * k)*(epsGamma * p) / (k * p);
            vertex(iGamma, iPsi, amp);
        };
    };
}

void EvtSVP_::decay_3body(EvtParticle* root) {
    root ->initializePhaseSpace(getNDaug(), getDaugs());

    EvtVector4R p = root->getDaug(0)->getP4(), // J/psi momentum
            k1 = root->getDaug(1)->getP4(), // mu+ momentum
            k2 = root->getDaug(2)->getP4(), // mu- momentum
            k = k1 + k2; // photon momentum

    double kSq = k*k;
    if (kSq < 1e-10) {
        return;
    }
    double kp = k*p;
    if (fabs(kp) < 1e-10) {
        return;
    }

    double dSq = delta*delta;
    double dSqDenom = dSq - k.mass2();
    if (fabs(dSqDenom) < 1e-10) {
        return;
    }

    double factor = dSq / (dSqDenom * kSq);

    for (int iPsi = 0; iPsi < 3; ++iPsi) {
        EvtVector4C epsPsi = root->getDaug(0)->epsParent(iPsi).conj();
        for (int iMplus = 0; iMplus < 2; ++iMplus) {
            EvtDiracSpinor spMplus = root->getDaug(1)->spParent(iMplus);
            for (int iMminus = 0; iMminus < 2; ++iMminus) {
                EvtDiracSpinor spMminus = root->getDaug(2)->spParent(iMminus);
                EvtVector4C epsGamma = EvtLeptonVCurrent(spMplus, spMminus);
                EvtComplex amp = (epsPsi * epsGamma) - (epsPsi * k)*(epsGamma * p) / kp;
                amp *= factor;
                vertex(iPsi, iMplus, iMminus, amp);
            };
        };
    };
}

void EvtSVP_::decay(EvtParticle *root) {
    if (getNDaug() == 2) decay_2body(root);
    else if (getNDaug() == 3) decay_3body(root);
}

void EvtSVP_::init() {
    checkSpinParent(EvtSpinType::SCALAR);
    if (getNDaug() == 2) { // old SVP
        checkNArg(0);
        checkNDaug(2);
        checkSpinDaughter(0, EvtSpinType::PHOTON);
        checkSpinDaughter(1, EvtSpinType::VECTOR);
    } else if (getNDaug() == 3) { // chi -> psi mu mu
        checkSpinParent(EvtSpinType::SCALAR);
        checkSpinDaughter(0, EvtSpinType::VECTOR);
        checkSpinDaughter(1, EvtSpinType::DIRAC);
        checkSpinDaughter(2, EvtSpinType::DIRAC);
        checkNArg(1);
        delta = getArg(0);
    }

}

void EvtSVP_::initProbMax() {
    if (getNDaug() == 2) setProbMax(2.5);
    else if (getNDaug() == 3) {
        if (getDaug(1).getId() == EvtPDL::getId("mu+").getId()) {
            setProbMax(550); // tested on 1e6 events
        }
        if (getDaug(1).getId() == EvtPDL::getId("e+").getId()) {
            setProbMax(8e3); // tested on 1e6 events
        }
    };
}
