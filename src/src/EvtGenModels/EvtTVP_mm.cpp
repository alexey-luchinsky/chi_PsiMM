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
// Module: EvtTVP_mm.cc
//
// Description: Routine to implement radiative decay chi_c2 -> psi gamma
//			matrix element from [S.P Baranov et al, PRD 85, 014034 (2012)]
//
// Modification history:
//	AVL	6 July, 2012	Module created
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


#include "EvtGenModels/EvtTVP_mm.hh"


#include <string>
#include <iostream>

using namespace std;



EvtTVP_mm::~EvtTVP_mm() {
}

EvtDecayBase* EvtTVP_mm::clone(){
  return new EvtTVP_mm;
}

void EvtTVP_mm::decay( EvtParticle *root ){
  root ->initializePhaseSpace(getNDaug(),getDaugs());

  EvtVector4R p=root->getDaug(0)->getP4(), // J/psi momentum
    k1 = root->getDaug(1)->getP4(),        // mu+ momentum
    k2 = root->getDaug(2)->getP4(),        // mu- momentum
    k=k1+k2;                               // photon momentum
  int iPols[4];
  for(int iChi= 0; iChi<5; iChi++) {
    iPols[0]=iChi;
    EvtTensor4C epsChi = root->epsTensor(iChi);
    for(int iPsi = 0; iPsi < 3; iPsi++) {
      iPols[1]=iPsi;
      EvtVector4C epsPsi = root->getDaug(0)->epsParent(iPsi).conj();
      for(int iMplus=0; iMplus<2; ++iMplus) {
        iPols[2]=iMplus;
        EvtDiracSpinor spMplus=root->getDaug(1)->spParent(iMplus).conj();
        for(int iMminus=0; iMminus<2; ++iMminus) {
          iPols[3]=iMminus;
          EvtDiracSpinor spMminus=root->getDaug(2)->spParent(iMminus);
          EvtVector4C epsGamma=EvtLeptonVCurrent(spMplus,spMminus);

          // [Baranov, (11)
          // matr = p^mu epsPsi^a epsChi_{a b} ( k_mu epsGamma_b  - k_b epsGamma_mu

          EvtVector4C eee = epsChi.cont1(epsPsi);
          EvtVector4C vvv = (p*k)*eee - (k*eee)*p;
          EvtComplex amp = vvv*epsGamma;
          amp = amp/(k*k);
          if(k.mass2()<0.0005) amp=0;
          vertex(iPols, amp);
        };
      };
    };
  };
}

void EvtTVP_mm::init(){
  checkNArg(0);
  checkNDaug(3);
  checkSpinParent(EvtSpinType::TENSOR);
  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::DIRAC);
}

void EvtTVP_mm::initProbMax() {
  if(getDaug(1).getId() == EvtPDL::getId("mu+").getId()) {
    cout<<"mu+"<<endl;
    setProbMax(520);  // tested on 1e6 events
  }
  if(getDaug(1).getId() == EvtPDL::getId("e+").getId()) {
    cout<<"e+"<<endl;
    setProbMax(2e6); // tested on 1e6 events
  }
  else {
    cout<<" EvtID "<<getDaug(1)<<" not realized yet"<<endl;
  }
}
