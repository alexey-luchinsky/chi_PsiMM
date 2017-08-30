#include <iostream>
#include <iomanip>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "EvtGenBase/EvtVector4R.hh"
#include "TLorentzVector.h"
#include "TH1F.h"

using namespace std;

void write_histogram_to_file(TH1F &histogram, string file_name) {
    const char *__file_name__ = file_name.c_str();
    remove(__file_name__);
    ofstream file;
        cout<<"Extporting histogram "<<file_name<<endl;
    file.open(__file_name__);
    for (int i = 1; i <= histogram.GetNbinsX(); i++)
        file << setiosflags(ios::scientific) << histogram.GetBinCenter(i) <<
        " " << setiosflags(ios::scientific) << histogram.GetBinContent(i) / histogram.GetBinWidth(i) <<
        " " << setiosflags(ios::scientific) << histogram.GetBinError(i) / histogram.GetBinWidth(i) << endl;

    file.close();
}


int main(int argc, char const *argv[]) {
  const char *root_file_name;
  if(argc>1)
    root_file_name=argv[1];
  else {
    cout<<"format: ./read_data.exe inFileName"<<endl;
    return 1;
  }
  TFile file(root_file_name);
  TTree *moms=(TTree*)file.Get("moms");
  cout<<moms->GetEntries()<<" entries in the tree"<<endl;
  TLorentzVector *_pPsi, *_k1, *_k2, *_kk1, *_kk2;
  double q2;
  moms->SetBranchAddress("q2",&q2);
  moms->SetBranchAddress("pPsi",&_pPsi);
  moms->SetBranchAddress("k1",&_k1);
  moms->SetBranchAddress("k2",&_k2);
  moms->SetBranchAddress("kk1",&_kk1);
  moms->SetBranchAddress("kk2",&_kk2);

  double minQ2=moms->GetMinimum("q2"), maxQ2=moms->GetMaximum("q2");
  int nBins=20;
  TH1F *hQ2=new TH1F("hQ2","hQ2",nBins,minQ2,maxQ2);
  for(int iEv=0; iEv<moms->GetEntries(); ++iEv) {
    moms->GetEntry(iEv);
    hQ2->Fill(q2);
  };
  write_histogram_to_file(*hQ2,"hh.hst");
  return 0;
}
