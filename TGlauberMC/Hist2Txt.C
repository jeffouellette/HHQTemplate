#include "runglauber_v3.1.C"

const int nbins = 200; //TODO: change me
const double max_x = 5; // TODO: change me

const bool saveQuarkDists = true; // TODO: whether to save wounded nucleon distributions in individual root files.

// Output file from runglauber should be called "outFile_${system}_${firstEvent}_${lastEvent}.root".
const string system = "pPb"; //TODO: e.g. "He3Au" or "pPb". For accessing and saving the right files/directories.
const int firstEvent = 500;
const int lastEvent = 750;

// TODO: select an energy scaling value
const double e0 = 0.00150022 * TMath::Power(140.*max_x/(5.*nbins), 4); // value for minbias pPb at 8TeV, smearing of 0.4 (better)

//const double e0_init = 0.0009;
//const double de0 = 0.0001;

void Hist2Txt () {
  // inFile contains the histogram with the wounded nucleon distribution for all events.
  TFile* inFile = new TFile(Form("outFile_%s_%i_%i.root", system.c_str(), firstEvent, lastEvent), "UPDATE");
  TFile* outRootFile = NULL;

  double max = 0; // for finding the maximum 
  for (int eventNum = firstEvent; eventNum < lastEvent; eventNum++) {
    //e0 = (e0_init) + de0 * (int)((eventNum - firstEvent) / 100);
    //cout << "event " << eventNum << ", e0=" << e0 << endl;

    TH2D* initedHist = (TH2D*)inFile->Get(Form("inited_event%i_translated", eventNum));
    initedHist->Scale (e0);

    ofstream outFile;
    char* outFileName = Form("./initedFiles_%s/event.dat", system.c_str());
    sprintf(outFileName, Form("./initedFiles_%s/event%d.dat", system.c_str(), eventNum));
    outFile.open(outFileName);

    if (saveQuarkDists) {
     // outRootFile copies the histogram with the distribution of wounded nucleons for a single event.
     // This can be used by diffusion.C in generating the initial positions of heavy qqbar pairs.
     TFile* outRootFile = new TFile(Form("./initedFiles_%s/event%d.root", system.c_str(), eventNum), "RECREATE");
     outRootFile->cd();
     initedHist->Write();
    }

    for (int xbin = 1; xbin <= nbins; xbin++) {
      for (int ybin = 1; ybin <= nbins; ybin++) {
        const long double content = initedHist->GetBinContent(xbin, ybin); 

        if (content == 0.) // if a value = 0 somewhere then hydro will fail!
          cout << "Found a zero in event " << eventNum << "!!!" << endl; // if this line prints at all you need to recreate the initial condition without 0's.

        outFile << content << "\t";
        if (content > max) max = content;
      }
      outFile << "\n";
    }
    outFile.close();
  }
  cout << "Max val: " << max << " GeV" << endl; // prints out the maximum initial energy density found

  inFile->Close();
  delete inFile;
}
