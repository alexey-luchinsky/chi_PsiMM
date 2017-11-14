#include "EvtGen/EvtGen.hh"
#include <string.h>

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

using namespace std;

void setLorentzVector(TLorentzVector *out, EvtVector4R in) {
    out->SetXYZT(in.get(1), in.get(2), in.get(3), in.get(0));
}

string to_str(double dbl) {
    char buffer[100];
    sprintf(buffer, "%f", dbl);
    return buffer;
};

int main(int argc, char** argv) {
    // ======== READ params ==========================
    cout << "Format: ./chic_all.exe inParticle [nEv=1e6] [decay_file=my_decay.dec] [postfix=]" << endl;
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


    string outName = string("root_") + string(argv[1]) + postfix + string(".root");
    TFile file(outName.c_str(), "RECREATE");
    TNtuple tup("tup", "tup", "id:q2:m2PsiK1:m2PsiK2:cosThEE:Mchi:m2K1KK1");
    TTree *moms = new TTree("moms", "moms");
    EvtVector4R pPsi, k1, k2, kk1, kk2;
    TLorentzVector *_pPsi = new TLorentzVector;
    moms->Branch("pPsi", "TLorentzVector", &_pPsi);
    TLorentzVector *_k1 = new TLorentzVector;
    moms->Branch("k1", "TLorentzVector", &_k1);
    TLorentzVector *_k2 = new TLorentzVector;
    moms->Branch("k2", "TLorentzVector", &_k2);
    TLorentzVector *_kk1 = new TLorentzVector;
    moms->Branch("kk1", "TLorentzVector", &_kk1);
    TLorentzVector *_kk2 = new TLorentzVector;
    moms->Branch("kk2", "TLorentzVector", &_kk2);
    double Q2;
    moms->Branch("q2", &Q2);

    // ======== MAIN LOOP ==========================
    int i;
    for (i = 0; i < nEvents; i++) {
        if (i % (nEvents / 10) == 0) {
            cout << "========= " << (int) (100. * i / nEvents) << " % ==========" << endl;
        };
        // Set up the parent particle
        EvtVector4R pInit(EvtPDL::getMass(CHI), 0.0, 0.0, 0.0);
        parent = EvtParticleFactory::particleFactory(CHI, pInit);
        parent->setDiagonalSpinDensity();

        // Generate the event
        myGenerator->generateDecay(parent);
        //        cout << "[DEBUG] chic_all: 1" << endl;
        // save the results
        EvtParticle *psi = parent->getDaug(0);
        pPsi = psi->getP4Lab();
        setLorentzVector(_pPsi, pPsi);
        if (i < 10)
            cout << " i=" << i << " pPsi=" << pPsi << endl;
        double m2PsiK1, m2PsiK2, cosThEE;
        EvtVector4R Ptot;
        k1 = parent->getDaug(1)->getP4Lab();
        setLorentzVector(_k1, k1);
        k2 = parent->getDaug(2)->getP4Lab();
        setLorentzVector(_k2, k2);
        EvtVector4R k = k1 + k2;
        if (psi->getNDaug() == 3) {
            kk1 = psi->getDaug(0)->getP4Lab();
            setLorentzVector(_kk1, kk1);
            kk2 = psi->getDaug(1)->getP4Lab();
            setLorentzVector(_kk2, kk2);
        };

        Ptot = pPsi + k1 + k2;
        Q2 = k*k;
        m2PsiK1 = (pPsi + k1)*(pPsi + k1);
        m2PsiK2 = (pPsi + k2)*(pPsi + k2);
        cosThEE = k.get(3) / k.d3mag();
        tup.Fill(parent->getPDGId(), Q2, m2PsiK1, m2PsiK2, cosThEE, sqrt(Ptot * Ptot), (k1 + kk1).mass2());
        moms->Fill();
        if (i < 0) {
            cout << "(* debug print at i=" << i << "========== *)" << endl;
            cout << " pPsi=" << pPsi << endl;
            cout << " Mpsi=" << sqrt(pPsi * pPsi) << endl;
            cout << " k1=" << k1 << endl;
            cout << " mmu=" << sqrt(k1 * k1) << endl;
            cout << " k2=" << k2 << endl;
            cout << " mmu=" << sqrt(k2 * k2) << endl;
            cout << " Ptot=" << Ptot << endl;
            cout << " Q2=" << (k1 + k2)*(k1 + k2) << endl;
            cout << " cosThEE=" << cosThEE << endl;
            cout << "  psi daugs=" << parent->getDaug(0)->getNDaug() << endl;
            cout << " kk1=" << parent->getDaug(0)->getDaug(0)->getP4Lab().mass() << endl;
            cout << " Mtot=" << (kk1 + kk2 + k1 + k2).mass() << endl;
        };
        parent->deleteTree();
    }

    tup.Write();
    moms->Write();
    file.Save();

    return 0;
}
