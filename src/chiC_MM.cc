#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"

#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtSimpleRandomEngine.hh"
#include "EvtGenBase/EvtMTRandomEngine.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtModel.hh"

#include "EvtGenModels/EvtSVP_mm.hh"
#include "EvtGenModels/EvtVVP_mm.hh"
#include "EvtGenModels/EvtTVP_mm.hh"

EvtGen *myGenerator;


void init_evtgen()
{
  EvtParticle* parent(0);
  EvtRandomEngine* eng = 0;
  #ifdef EVTGEN_CPP11
  // Use the Mersenne-Twister generator (C++11 only)
  eng = new EvtMTRandomEngine();
#else
  eng = new EvtSimpleRandomEngine();
#endif

  EvtRandom::setRandomEngine(eng);


  EvtAbsRadCorr* radCorrEngine = 0;

  std::list<EvtDecayBase*> extraModels;
  EvtModel &modelist=EvtModel::instance();
  modelist.registerModel(new EvtSVP_mm);
  modelist.registerModel(new EvtVVP_mm);
  modelist.registerModel(new EvtTVP_mm);

  myGenerator=new EvtGen("../src/my.dec" ,"../src/evt.pdl", eng,
                     radCorrEngine, &extraModels);


}



void write_histogram_to_file(TH1F &histogram, string file_name) {
    const char *__file_name__ = file_name.c_str();
    remove(__file_name__);
    ofstream file;

    file.open(__file_name__);
    for (int i = 1; i <= histogram.GetNbinsX(); i++)
        file << setiosflags(ios::scientific) << histogram.GetBinCenter(i) <<
        " " << setiosflags(ios::scientific) << histogram.GetBinContent(i) / histogram.GetBinWidth(i) <<
        " " << setiosflags(ios::scientific) << histogram.GetBinError(i) / histogram.GetBinWidth(i) << endl;

    file.close();
}

int main(int argc, char const *argv[])
{
  cout<<"Format: ./chic_all.exe inParticle [nEv=1e6] [decay_file=my_decay.dec]"<<endl;
   if(argc<2) {
     cout<<"Wrong number of arguments!"<<endl;
     return 1;
   };
   cout<<" Decaying particle "<<argv[1]<<endl;

   int nEvents(1e6);
   if(argc>2) nEvents=(int)atof(argv[2]);
   cout<<" nEvents="<<nEvents<<endl;



   char *decay_file;
//   if(argc>3) decay_file=argv[3];
//   else decay_file=(char*)"my_decay.dec";
  decay_file=(char*)"my_decay.dec";
   cout<<" decay_file = "<<decay_file<<endl;

  init_evtgen();

  EvtParticle* parent(0);

  static EvtId CHI = EvtPDL::getId(std::string(argv[1]));


  string outName=string("root_")+string(argv[1])+string(".root");
  TFile file(outName.c_str(),"RECREATE");
  TNtuple tup("tup","tup","id:q2:m2PsiK1:m2PsiK2:cosThEE:Mchi");
  // Loop to create nEvents, starting from an Upsilon(4S)
  int i;
  for (i = 0; i < nEvents; i++) {
    if(i%(nEvents/10)==0) {
      cout<<"========= "<<(int)(100.*i/nEvents)<<" % =========="<<endl;
    };
    // Set up the parent particle
    EvtVector4R pInit(EvtPDL::getMass(CHI), 0.0, 0.0, 0.0);
    parent = EvtParticleFactory::particleFactory(CHI, pInit);
    parent->setDiagonalSpinDensity();

    // Generate the event
    myGenerator->generateDecay(parent);
    // save the results
    EvtVector4R pPsi = parent->getDaug(0)->getP4Lab();
    EvtVector4R k1=parent->getDaug(1)->getP4Lab();
    EvtVector4R k2=parent->getDaug(2)->getP4Lab();
    EvtVector4R k = k1+k2;
    EvtVector4R Ptot=pPsi+k1+k2;
    double Q2=k*k;
    double m2PsiK1 = (pPsi+k1)*(pPsi+k1);
    double m2PsiK2 = (pPsi+k2)*(pPsi+k2);
    double cosThEE=k.get(3)/k.d3mag();
    tup.Fill(parent->getPDGId(), Q2, m2PsiK1, m2PsiK2, cosThEE, sqrt(Ptot*Ptot));
    if(i<0) {
      cout<<"(* debug print at i="<<i<<"========== *)"<<endl;
      cout<<" pPsi="<<pPsi<<endl;
      cout<<" Mpsi="<<sqrt(pPsi*pPsi)<<endl;
      cout<<" k1="<<k1<<endl;
      cout<<" mmu="<<sqrt(k1*k1)<<endl;
      cout<<" k2="<<k2<<endl;
      cout<<" mmu="<<sqrt(k2*k2)<<endl;
      cout<<" Ptot="<<Ptot<<endl;
      cout<<" Q2="<<(k1+k2)*(k1+k2)<<endl;
      cout<<" cosThEE="<<cosThEE<<endl;
    };
    parent->deleteTree();

  }

//  delete eng;
  tup.Write();
  file.Save();

  return 0;
}
