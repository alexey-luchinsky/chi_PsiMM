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
// Module: EvtSVP_mm.cc
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
#include "EvtGenBase/EvtDiracSpinor.hh"

#include "EvtGenModels/EvtSVP_mm.hh"


#include <string>
#include <iostream>

using namespace std;



EvtSVP_mm::~EvtSVP_mm() {
}

std::string EvtSVP_mm::getName(){
  return "SVP_mm";
}


EvtDecayBase* EvtSVP_mm::clone(){
  return new EvtSVP_mm;

}

void EvtSVP_mm::decay( EvtParticle *root ){
  root ->initializePhaseSpace(getNDaug(),getDaugs());

  EvtVector4R p=root->getDaug(0)->getP4(), // J/psi momentum
    k1 = root->getDaug(1)->getP4(),        // mu+ momentum
    k2 = root->getDaug(2)->getP4(),        // mu- momentum
    k=k1+k2;                               // photon momentum

  double kSq = k*k;
  if (kSq < 1e-10) {return;}
  double kp = k*p;
  if (fabs(kp) < 1e-10) {return;}

  double dSq = delta*delta;
  double dSqDenom = dSq - k.mass2();
  if (fabs(dSqDenom) < 1e-10) {return;}

  double factor = dSq/(dSqDenom*kSq);

  for(int iPsi=0; iPsi<3; ++iPsi) {
    EvtVector4C epsPsi = root->getDaug(0)->epsParent(iPsi).conj();
    for(int iMplus=0; iMplus<2; ++iMplus) {
      EvtDiracSpinor spMplus=root->getDaug(1)->spParent(iMplus);
      for(int iMminus=0; iMminus<2; ++iMminus) {
        EvtDiracSpinor spMminus=root->getDaug(2)->spParent(iMminus);
        EvtVector4C epsGamma=EvtLeptonVCurrent(spMplus,spMminus);
        EvtComplex amp = (epsPsi*epsGamma) - (epsPsi*k)*(epsGamma*p)/kp;
        amp *= factor;
        vertex(iPsi, iMplus, iMminus, amp);
      };
    };
  };
}


void EvtSVP_mm::init(){
  checkSpinParent(EvtSpinType::SCALAR);
  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::DIRAC);
  checkNArg(1);
  delta = getArg(0);
  checkNDaug(3);
  cout<<"EvtSVP::init() delta="<<delta<<endl;
}

void EvtSVP_mm::initProbMax() {
    if(getDaug(1).getId() == EvtPDL::getId("mu+").getId()) {
      setProbMax(550); // tested on 1e6 events
    }
    if(getDaug(1).getId() == EvtPDL::getId("e+").getId()) {
      setProbMax(8e3); // tested on 1e6 events
    }
    else {
      cout<<" EvtID "<<getDaug(1)<<" not realized yet"<<endl;
    }
}
