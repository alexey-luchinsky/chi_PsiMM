#ifndef EVTVVP_mm_HH
#define EVTVVP_mm_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtVVP_mm:public  EvtDecayAmp  {

public:

  EvtVVP_mm() {}
  virtual ~EvtVVP_mm();
  
  std::string getName();
  EvtDecayBase* clone();

  void initProbMax();
  void init();
  void decay(EvtParticle *p);
private:
  double delta;

};

#endif
