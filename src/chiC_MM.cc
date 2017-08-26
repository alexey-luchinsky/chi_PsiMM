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


void init_evtgen(int n=3)
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
  cout<<"Format: evtgen.exe [N=3]"<<endl;
  int N=3;
  if(argc>1) N=atoi(argv[1]);
  if(N !=3 && N!=5) {
    cout<<"N="<<N<<" is not realized yet"<<endl;
    return -1;
  }
  cout<<"N="<<N<<endl;
  TFile file(("../src/evtgen/omegaQQ_"+to_string(N)+"pi.root").c_str(),"RECREATE");
  TNtuple tup("tup","tup","Q2");

  init_evtgen(N);
  static EvtId CHI = EvtPDL::getId("Omega_bb+");
  EvtParticle* parent(0);
  int nEv=1e4;

  for(int iEv=0; iEv<nEv; ++iEv)
  {
    if( iEv % (nEv/10) == 0) cout<<"iEv="<<iEv<<endl;
    // Set up the parent particle
    EvtVector4R pInit(EvtPDL::getMass(CHI), 0.0, 0.0, 0.0);
    parent = EvtParticleFactory::particleFactory(CHI, pInit);
    parent->setDiagonalSpinDensity();
    // Generate the event
    myGenerator->generateDecay(parent);
    EvtVector4R pOout = parent->getDaug(0)->getP4Lab();
    double Q2=(pInit-pOout).mass2();
    tup.Fill(Q2);


    if(iEv<0)
    {
      cout<<pOout<<endl;
      cout<<"MOout="<<pOout.mass()<<endl;
      cout<<"==============="<<endl;
    }




  };
  tup.Write();
  TH1F hQ2("hQ2","hQ2",50,tup.GetMinimum("Q2"),tup.GetMaximum("Q2"));
  hQ2.Sumw2();
  tup.Project("hQ2","Q2");
  write_histogram_to_file(hQ2, "../src/evtgen/hQ2_"+to_string(N)+".txt");
  file.Save();




  return 0;
}
