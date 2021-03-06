TFile *_file;

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


void saveData(int id, string postfix, int ngroup=1) {
  string query;
  string name;
  switch(id) {
  case 10441:     query="id==10441";     name="chic0"; break;
  case 20443:     query="id==20443";     name="chic1"; break;
  case 445:      query="id==445";        name="chic2"; break;
  case 10551:     query="id==10551";     name="chib0"; break;
  case 20553:     query="id==20553";     name="chib1"; break;
  case 555:      query="id==555";        name="chib2"; break;
  default:      cout<<"Unknown ID"<<endl;    return;
  };
  name = name + postfix;
  cout<<"Exporting histograms for id="<<id<<" ("<<name<<"), "<<tup->GetEntries(query.c_str())<<" entries"<<endl;
  cout<<" query="<<query<<endl;
  if(tup->GetEntries(query.c_str())>0) {
    tup->Project("hMchi0","Mchi",query.c_str());
    write_histogram_to_file(*hMchi0,("hst/hM_"+name+".hst").c_str(),ngroup);
    hMchi0->Delete();

    tup->Project("hQ20","q2",query.c_str());
    write_histogram_to_file(*hQ20,("hst/hQ2_"+name+".hst").c_str(),ngroup);
    hQ20->Delete();

    tup->Project("hCos1","cosThEE",query.c_str());
    write_histogram_to_file(*hCos1,("hst/hCos_"+name+".hst").c_str(),ngroup);
    hCos1->Delete();

    tup->Project("hm2PsiK1","m2PsiK1",query.c_str());
    write_histogram_to_file(*hm2PsiK1,("hst/hm2PsiK1_"+name+".hst").c_str(),ngroup);
    hm2PsiK1->Delete();

    tup->Project("hm2K1KK1","m2K1KK1",query.c_str());
    write_histogram_to_file(*hm2K1KK1,("hst/hm2K1KK1_"+name+".hst").c_str(),ngroup);
    hm2K1KK1->Delete();
  }
  else {
    cout<<" No events for id="<<id<<endl;
  };
  cout<<"==========="<<endl;

}

TChain *chainc, *chainb;

read_data_chi(string fileName="max", string postfix="",int ngroup=1) {
  TChain chain("tup");
  chain.Add(fileName.c_str());
  cout<<tup->GetEntries()<<" entries"<<endl;
  int idChi0=10441, idChi1=20443, idChi2=445;
  int idChiB0=10551, idChiB1=20553, idChiB2=555;
  saveData(idChi0,postfix,ngroup);  saveData(idChi1,postfix,ngroup); saveData(idChi2,postfix,ngroup);
  saveData(idChiB0,postfix,ngroup); saveData(idChiB1,postfix,ngroup);saveData(idChiB2,postfix,ngroup);
}
