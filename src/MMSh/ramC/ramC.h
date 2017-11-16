#ifndef RAMC
#define RAMC

//#include "CLHEP/Vector/LorentzVector.h"
//#include "EvtGenBase/EvtVector4R.hh"
//#include "CLHEP/Random/Random.h"
//using namespace CLHEP;

class RAMBOC {
public:
	RAMBOC(int n, double ecm, double xm[]);
	double next();
	double getWT() {return WT;};
	double getECM() { return ECM;};
	void setECM(double ecm) { ECM=ecm;};
	double getMass(int i) { return XM[i];};
	int getN() {return N;};
	double* getP(int i);

private:
	double N, ECM, XM[5], WT;
	double p1[4], p2[4], p3[4], p4[4], p5[4];
};

#endif
