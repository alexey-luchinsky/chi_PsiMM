#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <iostream>
#include <string>
#include "EvtGenBase/EvtVector3C.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenModels/EvtVVP_mm.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtTensor4C.hh"

using namespace std;

EvtVVP_mm::~EvtVVP_mm() {}

std::string EvtVVP_mm::getName(){
  return "VVP_mm";
}


EvtDecayBase* EvtVVP_mm::clone(){
  return new EvtVVP_mm;
}

void EvtVVP_mm::init(){
  checkNDaug(3);
  checkSpinParent(EvtSpinType::VECTOR);
  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::DIRAC);
}


void EvtVVP_mm::decay( EvtParticle *root ){
  root ->initializePhaseSpace(getNDaug(),getDaugs());
  EvtVector4R
    pChi = root->getP4Lab(),                 // chi momentum
    pPsi = root->getDaug(0)->getP4Lab(),   // psi momentum
    k1 = root->getDaug(1)->getP4(),        // mu+ momentum
    k2 = root->getDaug(2)->getP4(),        // mu- momentum
    k=k1+k2;                               // photon momentum
  int iPols[4];
  for(int iChi= 0; iChi<3; iChi++) {
    iPols[0]=iChi;
    EvtVector4C epsChi = root->epsParent(iChi);
    for(int iPsi = 0; iPsi < 3; iPsi++) {
      iPols[1]=iPsi;
      EvtVector4C epsPsi = root->getDaug(0)->epsParent(iPsi).conj();
      for(int iMplus=0; iMplus<2; ++iMplus) {
        iPols[2]=iMplus;
        EvtDiracSpinor spMplus=root->getDaug(1)->spParent(iMplus);
        for(int iMminus=0; iMminus<2; ++iMminus) {
          iPols[3]=iMminus;
          EvtDiracSpinor spMminus=root->getDaug(2)->spParent(iMminus);
          EvtVector4C epsGamma=EvtLeptonVCurrent(spMplus,spMminus).conj();
          // amp = e_{mu nu alpha beta} epsChi^mu epsPsi^nu epsGamma^alpha k^beta/k^2
          EvtComplex amp = k*dual(EvtGenFunctions::directProd(epsChi,epsPsi)).cont1(epsGamma);
          amp = amp/(k*k);
          if(k.mass2()<0.0005) amp=0;
          vertex(iPols, amp);
        };
      };
    };
  };
}

void EvtVVP_mm::initProbMax(){
  if(getDaug(1).getId() == EvtPDL::getId("mu+").getId()) {
    setProbMax(65); // tested on 1e6 events
  }
  if(getDaug(1).getId() == EvtPDL::getId("e+").getId()) {
    setProbMax(1e3); // tested on 1e6 events
  }
  else {
    cout<<" EvtID "<<getDaug(1)<<" not realized yet"<<endl;
  }
}
