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


#include "EvtGenModels/EvtVVP_.hh"


#include <string>
#include <iostream>

using namespace std;

EvtVVP_::~EvtVVP_() {

}

std::string EvtVVP_::getName() {
    return "VVP_";
}

EvtDecayBase* EvtVVP_::clone() {
    return new EvtVVP_;

}

void EvtVVP_::decay_2body(EvtParticle* p) {

    p->initializePhaseSpace(getNDaug(), getDaugs());

    EvtParticle *v, *ph;

    v = p->getDaug(0);
    ph = p->getDaug(1);

    EvtVector3C epsp[3];
    EvtVector3C epsv[3];
    EvtVector3C epsph[2];

    epsp[0] = p->eps(0).vec();
    epsp[1] = p->eps(1).vec();
    epsp[2] = p->eps(2).vec();

    epsv[0] = v->eps(0).vec().conj();
    epsv[1] = v->eps(1).vec().conj();
    epsv[2] = v->eps(2).vec().conj();

    epsph[0] = ph->epsParentPhoton(0).vec().conj();
    epsph[1] = ph->epsParentPhoton(1).vec().conj();

    int i, j, k;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 2; k++) {
                vertex(i, j, k, epsp[i].cross(epsv[j]) * epsph[k]);
            }
        }
    }

    return;
}

void EvtVVP_::decay_3body(EvtParticle* root) {
    root ->initializePhaseSpace(getNDaug(), getDaugs());
    EvtVector4R
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
    for (int iChi = 0; iChi < 3; iChi++) {
        iPols[0] = iChi;
        EvtVector4C epsChi = root->epsParent(iChi);
        for (int iPsi = 0; iPsi < 3; iPsi++) {
            iPols[1] = iPsi;
            EvtVector4C epsPsi = root->getDaug(0)->epsParent(iPsi).conj();
            for (int iMplus = 0; iMplus < 2; ++iMplus) {
                iPols[2] = iMplus;
                EvtDiracSpinor spMplus = root->getDaug(1)->spParent(iMplus);
                for (int iMminus = 0; iMminus < 2; ++iMminus) {
                    iPols[3] = iMminus;
                    EvtDiracSpinor spMminus = root->getDaug(2)->spParent(iMminus);
                    EvtVector4C epsGamma = EvtLeptonVCurrent(spMplus, spMminus).conj();
                    // amp = e_{mu nu alpha beta} epsChi^mu epsPsi^nu epsGamma^alpha k^beta/k^2
                    EvtComplex amp = k * dual(EvtGenFunctions::directProd(epsChi, epsPsi)).cont1(epsGamma);
                    amp *= factor;
                    vertex(iPols, amp);
                };
            };
        };
    };
}

void EvtVVP_::decay(EvtParticle *root) {
    if (getNDaug() == 2) decay_2body(root);
    else if (getNDaug() == 3) decay_3body(root);
}

void EvtVVP_::init() {
    checkSpinParent(EvtSpinType::VECTOR);
    if (getNDaug() == 2) { // old SVP
        checkNArg(0);
        checkNDaug(2);
        checkSpinDaughter(0, EvtSpinType::VECTOR);
        checkSpinDaughter(1, EvtSpinType::PHOTON);
    } else if (getNDaug() == 3) { // chi -> psi mu mu
        checkSpinDaughter(0, EvtSpinType::VECTOR);
        checkSpinDaughter(1, EvtSpinType::DIRAC);
        checkSpinDaughter(2, EvtSpinType::DIRAC);
        checkNArg(1);
        delta = getArg(0);
    }

}

void EvtVVP_::initProbMax() {
    if (getNDaug() == 2) setProbMax(1.5);
    else if (getNDaug() == 3) {
        if (getDaug(1).getId() == EvtPDL::getId("mu+").getId()) {
            setProbMax(15); // tested on 1e6 events
        }
        if (getDaug(1).getId() == EvtPDL::getId("e+").getId()) {
            setProbMax(8e3); // tested on 1e6 events
        }
    };
}

