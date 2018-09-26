
void plot_dNch_deta () {
  ifstream infile;
  infile.open("dNch_deta.dat");
  string dummyLine;
  getline(infile, dummyLine);

  int event;
  double dNch_dy;
  double e0;
  TH2F* dN_dNch = new TH2F("dN_dNch", "", 5, 0.00155, 0.00205, 50, 0, 150);
  /*for (int i=0;i<10;i++) {
    dN_dNch[i] = new TH1F*(Form("dN_dNch_e0bin_%i", i), "", 35, 0, 70);
  }*/
  while(infile) {
    infile >> event;
    infile >> e0;
    infile >> dNch_dy;
    //int e0bin = (int)(event / 100); 
    dN_dNch->Fill(e0, dNch_dy);
  }
  TCanvas* c = new TCanvas("c", "", 800, 600);

  dN_dNch->ProfileX()->Draw("e1");
  TProfile* profx = dN_dNch->ProfileX();

  TF1* linfit = new TF1("linfit", "[0] + [1]*x", -0.5, 9.5);
  profx->Fit(linfit);

  const float b = linfit->GetParameter(0);
  const float m = linfit->GetParameter(1);

  cout << "Optimal x-axis value = " << (20.1 - b) / m << endl;
 

  c->SaveAs("dNch_deta.pdf");
}
