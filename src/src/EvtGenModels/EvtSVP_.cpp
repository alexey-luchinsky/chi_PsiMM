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
//	AVL	Jul 6, 2012	modle created
//
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

void EvtSVP_::decay(EvtParticle *root) {
//    cout << "[DEBUG] EvtSVP_::decay()" << endl;
    ncall++;
    root ->initializePhaseSpace(getNDaug(), getDaugs());

    EvtVector4R p = root->getDaug(1)->getP4(), // J/psi momentum
            k = root->getDaug(0)->getP4(); // Photon momentum
    for (int iPsi = 0; iPsi < 4; iPsi++) {
        for (int iGamma = 0; iGamma < 1; iGamma++) {
//            cout << "[DEBUG] EvtSVP_::decay(): 1" << endl;
            EvtVector4C epsPsi = root->getDaug(1)->epsParent(iPsi).conj();
//            cout << "[DEBUG] EvtSVP_::decay(): 2" << endl;
            EvtVector4C epsGamma = root->getDaug(0)->epsParentPhoton(iGamma).conj();
//            cout << "[DEBUG] EvtSVP_::decay(): 3" << endl;

            EvtComplex amp = (epsPsi * epsGamma) - (epsPsi * k)*(epsGamma * p) / (k * p);
//            cout << "[DEBUG] EvtSVP_::decay(): 4" << endl;


            vertex(iGamma, iPsi, amp);
//            cout << "[DEBUG] EvtSVP_::decay(): 5" << endl;
        }
    }

}

void EvtSVP_::init() {

    ncall = 0;

    checkNArg(0);
    checkNDaug(2);


    checkSpinParent(EvtSpinType::SCALAR);

    checkSpinDaughter(0, EvtSpinType::PHOTON);
    checkSpinDaughter(1, EvtSpinType::VECTOR);

}

void EvtSVP_::initProbMax() {
    setProbMax(1.2);
}
