#include "EvtGen/EvtGen.hh"
#include <string.h>
#include <iomanip>

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
#include "EvtGenModels/EvtSVP_.hh"
#include "EvtGenModels/EvtVVP_.hh"
#include "EvtGenModels/EvtTVP_.hh"

#include "TFile.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "TH1F.h"
using namespace std;

void print_mean(TNtuple *tup, string var, string postfix);

void write_histogram_to_file(TH1F &histogram, string file_name, int ngroup=1) {
    const char *__file_name__ = file_name.c_str();
    remove(__file_name__);
    ofstream file;
        cout<<"Extporting histogram "<<file_name<<endl;
    file.open(__file_name__);
    histogram.Rebin(ngroup);
    for (int i = 1; i <= histogram.GetNbinsX(); i++)
        file << setiosflags(ios::scientific) << histogram.GetBinCenter(i) <<
        " " << setiosflags(ios::scientific) << histogram.GetBinContent(i) / histogram.GetBinWidth(i) <<
        " " << setiosflags(ios::scientific) << histogram.GetBinError(i) / histogram.GetBinWidth(i) << endl;
    file.close();
    cout<<"\t integral="<<histogram.Integral()<<" nBins="<<histogram.GetNbinsX()<<endl;
}


int main(int argc, char** argv) {
    // ======== READ params ==========================
    cout << "Format: ./test.exe inParticle [nEv=1e6] [decay_file=my_decay.dec] [postfix=]" << endl;
    if (argc < 2) {
        cout << "Wrong number of arguments!" << endl;
        return 1;
    };
    cout << " Decaying particle " << argv[1] << endl;

    int nEvents(1e6);
    if (argc > 2) nEvents = (int) atof(argv[2]);
    cout << " nEvents=" << nEvents << endl;


    char *decay_file;
    if (argc > 3) decay_file = argv[3];
    else decay_file = (char*) "../src/my.dec";
    cout << " decay_file = " << decay_file << endl;

    char *postfix;
    if (argc > 4) postfix = argv[4];
    else postfix = (char*) "";
    cout << "postfix=" << postfix << endl;

    // ======== INIT ==========================
    EvtParticle * parent(0);
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
    EvtModel &modelist = EvtModel::instance();
    modelist.registerModel(new EvtSVP_mm);
    modelist.registerModel(new EvtVVP_mm);
    modelist.registerModel(new EvtTVP_mm);

    modelist.registerModel(new EvtSVP_);
    modelist.registerModel(new EvtVVP_);
    modelist.registerModel(new EvtTVP_);

    EvtGen *myGenerator = new EvtGen(decay_file, "../src/evt.pdl", eng,
            radCorrEngine, &extraModels);



    static EvtId CHI = EvtPDL::getId(std::string(argv[1]));

    double e1,pz1,e2,pz2,e3,pz3;
    EvtVector4R p1,p2,p3;
    string outName = string("test_") + string(argv[1]) + postfix + string(".root");
//    TFile file(outName.c_str(), "RECREATE");
    TNtuple tup("tup", "tup", "e1:pz1:e2:pz2:e3:pz3");

    for (int iEv = 0; iEv < nEvents; iEv++) {
        if (iEv % (nEvents / 10) == 0) {
            cout << "========= " << (int) (100. * iEv / nEvents) << " % ==========" << endl;
        };
        // Set up the parent particle
        EvtVector4R pInit(EvtPDL::getMass(CHI), 0.0, 0.0, 0.0);
        parent = EvtParticleFactory::particleFactory(CHI, pInit);
        parent->setDiagonalSpinDensity();
        // Generate the event
        myGenerator->generateDecay(parent);
        p1=parent->getDaug(0)->getP4Lab();
        p2=parent->getDaug(1)->getP4Lab();
        if(parent->getNDaug()>2)
            p3=parent->getDaug(2)->getP4Lab();
        tup.Fill(p1.get(0),abs(p1.get(3)),p2.get(0),abs(p2.get(3)),p3.get(0),abs(p3.get(3)));
        parent->deleteTree();
    }

    tup.Write();
//    file.Save();
    cout<<"================"<<endl;
    string postfix_=string(argv[1])+"_"+postfix;
    print_mean(&tup,"e1",postfix_);print_mean(&tup,"pz1",postfix_);
    print_mean(&tup,"e2",postfix_);print_mean(&tup,"pz2",postfix_);
    print_mean(&tup,"e3",postfix_);print_mean(&tup,"pz3",postfix_);
    return 0;
}

void print_mean(TNtuple *tup, string var, string postfix) {
    const char *cvar=var.c_str();
    TH1F *hist=new TH1F(cvar,cvar,30,tup->GetMinimum(cvar),tup->GetMaximum(cvar));
    tup->Project(cvar,cvar);
    cout<<"<"<<var<<">="<<hist->GetMean()<<" "<<hist->GetRMS()<<endl;
    write_histogram_to_file(*hist,var+"_"+postfix+".hst");
    delete hist;
}
