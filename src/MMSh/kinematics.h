/* 
 * File:   kinematics.h
 * Author: luchinsky
 *
 * Created on November 16, 2017, 3:37 PM
 */

#ifndef KINEMATICS_H
#define	KINEMATICS_H
#include <iostream>
#include "math.h"

using namespace std;
void print_v4(double p[4]);
void print_v4(string name, double p[4]);
void println_v4(double p[4]);
void println_v4(string name, double p[4]);
double sp(double p1[4], double p2[4]);
double sp(double p[4]);
void sum(double *p1, double *p2, double *P);
void sum(double *p1, double *p2, double *p3, double *P);
double sum_mass2(double *p1, double *p2);
double sum_mass2(double *p1, double *p2, double *p3);
void apply_boost_to(double bx, double by, double bz, double (&answ)[4]);
void apply_boost_to(double (&boost)[4], double (&answ)[4]);


#endif	/* KINEMATICS_H */

