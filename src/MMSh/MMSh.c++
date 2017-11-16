#include <iostream>
#include "math.h"
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
    return ((double) rand() / (double) (RAND_MAX));
}

void print_v4(double p[4]) {
    cout << "{" << p[3] << "," << p[0] << "," << p[1] << "," << p[2] << "}";
}

void print_v4(string name, double p[4]) {
    cout << name << "=";
    print_v4(p);
}

void println_v4(double p[4]) {
    print_v4(p);
    cout << endl;
};

void println_v4(string name, double p[4]) {
    cout << name << "=";
    println_v4(p);
}

double sp(double p1[4], double p2[4]) {
    return p1[3] * p2[3] - p1[0] * p2[0] - p1[1] * p2[1] - p1[2] * p2[2];
}

double sp(double p[4]) {
    return sp(p, p);
}

void sum(double *p1, double *p2, double *P) {
    for (int i = 0; i < 4; ++i) P[i] = p1[i] + p2[i];
}

void sum(double *p1, double *p2, double *p3, double *P) {
    for (int i = 0; i < 4; ++i) P[i] = p1[i] + p2[i] + p3[i];
}

double sum_mass2(double *p1, double *p2) {
    double P[4];
    sum(p1, p2, P);
    return sp(P);
}

double sum_mass2(double *p1, double *p2, double *p3) {
    double P[4];
    sum(p1, p2, p3, P);
    return sp(P);
}
    double const PI = acos(-1), TWO_PI=2*PI;

void test_2body(int nEv) {
    const int nOut = 2;
    double mmu = 10;
    double XM[nOut] = {0, 0};
    double p1[4], p2[4];

    double sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = ram2_(mmu, XM, p1, p2);
        wt *= pow(TWO_PI,4-3*nOut);
        sum += wt;
    };
    sum /= nEv;
    cout<<" sum="<<sum<<" vs "<<1./(8*PI)<<endl;
    
}

void test_3body(int nEv) {
    const int nOut = 3;
    double mmu = 10;
    double XM[nOut] = {0, 0, 0};
    double p1[4], p2[4], p3[4], P[4];

    double sum = 0;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = ram3_(mmu, XM, p1, p2, p3);
        wt *= pow(TWO_PI,4-3*nOut);
        sum += wt;
    };
    sum /= nEv;
    cout<<" sum="<<sum<<endl;
}

int main(void) {
    test_2body(1e4);
}
