//
//  ramc_test.c++
//  
//
//  Created by Luchinsky on 20/10/14.
//
//

#include <stdio.h>
#include <TDirectoryFile.h>
#include "ramC.h"
#include "TFile.h"
#include "TNtuple.h"

extern "C" {
    double rndm2_(double dummy);
}

using namespace std;


void printVec(string name, HepLorentzVector v) {
    cout<<name<<"=(" <<v.getT()<<","<<v.getX()<<","<<v.getY()<<","<<v.getZ()<<");";
};

void printVecLn(string name, HepLorentzVector v) {
    printVec(name, v); cout<<endl;
};


int main(void) {
    TFile file("Test.root","RECREATE");
    TNtuple tup("tup","tup","s12:s13:s23:wt");
    
    srand (1);
    cout<<"Hello, world!"<<endl;
    const int N=3;
    double XM[N];
    XM[0]=0; XM[1]=0; XM[2]=0;
    RAMBOC ram(N,10,XM);
    for(int i=0; i<1000000; ++i) {
//        cout<<ram.getWT()<<endl;
        double wt=ram.next();
        HepLorentzVector p1=ram.getP(0),
                p2=ram.getP(1),
                p3=ram.getP(2);
        double s12=(p1+p2)*(p1+p2), s13=(p1+p3)*(p1+p3), s23=(p2+p3)*(p2+p3);
        tup.Fill(s12,s13,s23,wt);
        if(i<3) {
            cout<<"===="<<endl;
            printVecLn("p1",p1);
            cout<<"p1*p1="<<p1.m2()<<endl;
            printVecLn("p2",p2);
            cout<<"p2*p2="<<p2.m2()<<endl;
            printVecLn("p3",p3);
            cout<<"p3*p3="<<p3.m2()<<endl;
            printVecLn("Q",p1+p2+p3);
            cout<<"p1*p2="<<p1*p2<<endl;
        };
    };
    
    tup.Write(); file.Save(); file.Close();
    
    return 0;
}
