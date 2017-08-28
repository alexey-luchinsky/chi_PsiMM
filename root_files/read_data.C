TFile *_file;

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

void saveData(int id, string postfix) {
  string query;
  string name;
  if(id==10441) {
    query="id==10441";
    name="chic0";
  }
  else if(id==20443) {
    query="id==20443";
    name="chic1";
  }
  else if(id=445) {
    query="id==445";
    name="chic2";
  }
  else {
    cout<<"Unknown ID"<<endl;
    return;
  };
  cout<<"Exporting histograms for id="<<id<<" ("<<name<<"), "<<tup->GetEntries(query.c_str())<<" entries"<<endl;
  cout<<" query="<<query<<endl;
  if(tup->GetEntries(query.c_str())>0) {
    tup->Project("hMchi0","Mchi",query.c_str()); write_histogram_to_file(*hMchi0,("hM_"+name+".hst").c_str());
    tup->Project("hQ20","q2",query.c_str()); write_histogram_to_file(*hQ20,("hQ2_"+name+".hst").c_str());
    tup->Project("hCos1","cosThEE",query.c_str()); write_histogram_to_file(*hCos1,("hCos_"+name+".hst").c_str());
    tup->Project("hm2PsiK1","m2PsiK1",query.c_str()); write_histogram_to_file(*hm2PsiK1,("hm2PsiK1_"+name+".hst").c_str());
  }
  else {
    cout<<" No events for id="<<id<<endl;
  };
  cout<<"==========="<<endl;

}

read_data(char *file_name, string postfix="") {
        _file=new TFile(file_name);
        int idChi0=10441, idChi1=20443, idChi2=445;
      saveData(idChi0,postfix);
      saveData(idChi1,postfix);
      saveData(idChi2,postfix);
}
