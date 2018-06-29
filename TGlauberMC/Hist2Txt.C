#include "runglauber_v3.0.C"

const int nbins = 141; //TODO: change me
const string system = "He3Au"; //TODO: "He3Au" or "pPb"

void Hist2Txt(const int firstEvent = 0, // first event to process
              const int lastEvent  = 500) // last event to process
{
  // inFile contains the histogram with the wounded nucleon distribution for all events.
  TFile* inFile = new TFile(Form("outFile_%s.root", system.c_str()), "UPDATE");

  double max = 0; // for finding the maximum 
  for (int eventNum = firstEvent; eventNum < lastEvent; eventNum++) {
    TH2D* initedHist = (TH2D*)inFile->Get(Form("inited_event%i_translated", eventNum));
    ofstream outFile;
    char* outFileName = Form("./initedFiles_%s/event.dat", system.c_str());
    sprintf(outFileName, Form("./initedFiles_%s/event%d.dat", system.c_str(), eventNum));
    outFile.open(outFileName);

    // outRootFile copies the histogram with the distribution of wounded nucleons for a single event.
    // This can be used by diffusion.C in generating the initial positions of heavy qqbar pairs.
    TFile* outRootFile = new TFile(Form("./initedFiles_%s/event%d.root", system.c_str(), eventNum), "RECREATE");
    outRootFile->cd();
    initedHist->Write();

    for (int xbin = 1; xbin <= nbins; xbin++) { // loops exactly nbins+1 times (so NUMT in params should be EXACTLY nbins+1)
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
