#include "ramC.h"

#include <iostream>

using namespace std;

extern "C" {
	void init_(void);
	double ram2_(double &ecm, double *xm, double *p1, double *p2);
	double ram3_(double &ecm, double *xm, double *p1, double *p2, double *p3);
	double ram4_(double &ecm, double *xm, double *p1, double *p2, double *p3, double *p4);
	double ram5_(double &ecm, double *xm, double *p1, double *p2, double *p3, double *p4, double *p5);
    double rndm2_(double);
}



double rndm2_(double dummy) {
    return ((double)rand() / (double)(RAND_MAX));
}


RAMBOC::RAMBOC(int n, double ecm, double *xm) {
	N=n;
	ECM=ecm;
	for(int i=0; i<N; i++) XM[i]=xm[i];
        
//	init_();
};

double RAMBOC::next() {
	if(N==2) {
		WT=ram2_(ECM, XM, p1, p2);
	}
	else if(N==3) {
		WT=ram3_(ECM, XM, p1, p2, p3);
	}
	else if(N==4) {
		WT=ram4_(ECM, XM, p1, p2, p3,p4);
	}
	else if(N==5) {
		WT=ram5_(ECM, XM, p1, p2, p3,p4,p5);
	};
	return WT;
};

double*  RAMBOC::getP(int i) {
	if(i==0) return p1;
	else if(i==1) return p2;
	else if(i==2) return p3;
	else if(i==3) return p4;
	else if(i==4) return p5;
    else {
        cout<<"ramC: not supported yet"<<endl;
        abort();
    }
};

