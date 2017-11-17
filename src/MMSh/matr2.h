/* 
 * File:   matr2.h
 * Author: luchinsky
 *
 * Created on November 16, 2017, 5:29 PM
 */

#ifndef MATR2_H
#define	MATR2_H
extern double mmu, Mpsi;
extern double Mchi0, Mchi1, Mchi2;
extern double const PI;

double matr2_0(double pPsi[4], double k1[4], double k2[4]);
double _matr2LL_0(double kk1[4], double kk2[4], double k1[4], double k2[4]);
double matr2_1(double pPsi[4], double k1[4], double k2[4]);
double matr2_2(double pPsi[4], double k1[4], double k2[4]);

#endif	/* MATR2_H */

