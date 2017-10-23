
#ifndef EvtTVP_mm_HH
#define EvtTVP_mm_HH

#include <fstream>
#include <stdio.h>


#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtSemiLeptonicAmp.hh"

class EvtParticle;

class EvtTVP_mm:public  EvtDecayAmp  {

public:

  EvtTVP_mm(double _delta=1e5) {};
  virtual ~EvtTVP_mm();

  std::string getName() {return "TVP_mm";}
  EvtDecayBase* clone();

  void decay(EvtParticle *p);
  void init();

  virtual void initProbMax();
private:
  double delta;
};

#endif

