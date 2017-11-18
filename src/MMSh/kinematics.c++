#include "kinematics.h"

void print_v4(double p[4]) {
    cout << "{" << p[3] << "," << p[0] << "," << p[1] << "," << p[2] << "};";
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

void apply_boost_to(double bx, double by, double bz, double (&answ)[4]) {
    double bxx = bx * bx;
    double byy = by * by;
    double bzz = bz * bz;

    double b2 = bxx + byy + bzz;
    if (b2 == (double) 0.0)
        return;

    double gamma = 1.0 / sqrt(1 - b2);
    double gb2 = (gamma - 1.0) / b2;
    double gb2xy = gb2 * bx * by;
    double gb2xz = gb2 * bx * bz;
    double gb2yz = gb2 * by * bz;

    double gbx = gamma * bx;
    double gby = gamma * by;
    double gbz = gamma * bz;

    double e2 = answ[3];
    double px2 = answ[0];
    double py2 = answ[1];
    double pz2 = answ[2];

    answ[3] = gamma * e2 + gbx * px2 + gby * py2 + gbz * pz2;
    answ[0] = gbx * e2 + gb2 * bxx * px2 + px2 + gb2xy * py2 + gb2xz * pz2;
    answ[1] = gby * e2 + gb2 * byy * py2 + py2 + gb2xy * px2 + gb2yz * pz2;
    answ[2] = gbz * e2 + gb2 * bzz * pz2 + pz2 + gb2yz * py2 + gb2xz * px2;
}

void apply_boost_to(double (&boost)[4], double (&answ)[4]) {
    double bx = (double) (boost[0] / boost[3]);
    double by = (double) (boost[1] / boost[3]);
    double bz = (double) (boost[2] / boost[3]);
    apply_boost_to(bx, by, bz, answ);
}
