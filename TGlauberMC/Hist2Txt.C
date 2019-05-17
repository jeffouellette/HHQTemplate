#include <TMath.h>
#include <TFile.h>
#include <TH2.h>
#include <TNtuple.h>

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

const int nbins =  60; //TODO: change me
const double max_x = 6; // TODO: change me

const bool saveQuarkDists = true; // TODO: whether to save wounded nucleon distributions in individual root files.

// TODO: select an energy scaling value
const double e0 = 0.00150022 * TMath::Power (140.*max_x/(5.*nbins), 4); // value for minbias pPb at 8TeV, smearing of 0.4 (better)
//const double e0 = 0.375 * 0.00150022 * TMath::Power (140.*max_x/(5.*nbins), 4); // value for minbias dAu at 200GeV, smearing of 0.4

//const double e0_init = 0.0009;
//const double de0 = 0.0001;

void Hist2Txt (const int firstEvent, const int lastEvent, const char* sys1, const char* sys2) {
  // inFile contains the histogram with the wounded nucleon distribution for all events.
  TFile* inFile = new TFile (Form ("outputs/outFile_%s%s_%i_%i.root", sys1, sys2, firstEvent, lastEvent), "UPDATE");
  TFile* outRootFile = NULL;

  double max = 0; // for finding the maximum 
  char* outFileName;

  for (int eventNum = firstEvent; eventNum < lastEvent; eventNum++) {
    //e0 = (e0_init) + de0 * (int)((eventNum - firstEvent) / 100);
    //cout << "event " << eventNum << ", e0=" << e0 << endl;

    TH2D* h_ed = (TH2D*) inFile->Get (Form ("h_ed_event%i", eventNum));
    //h_ed = (TH2D*) h_ed->Clone (Form ("h_ed_event%i", eventNum));
    h_ed->Scale (e0);

    ofstream outFile;
    outFileName = Form ("./initedFiles_%s%s/event%d.dat", sys1, sys2, eventNum);
    outFile.open (outFileName);

    if (saveQuarkDists) {
      // outRootFile copies the histogram with the distribution of wounded nucleons for a single event.
      // This can be used by diffusion.C in generating the initial positions of heavy qqbar pairs.
      TFile* outRootFile = new TFile (Form ("./initedFiles_%s%s/event%d.root", sys1, sys2, eventNum), "RECREATE");
      outRootFile->cd ();
      h_ed->Write ();
    }

    for (int xbin = 1; xbin <= nbins; xbin++) {
      for (int ybin = 1; ybin <= nbins; ybin++) {
        const long double content = h_ed->GetBinContent (xbin, ybin); 

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

  TNtuple* nt = (TNtuple*) inFile->Get ("nt");
  outRootFile = new TFile (Form ("ntuples/tntuple_%s%s_%i_%i.root", sys1, sys2, firstEvent, lastEvent), "RECREATE");
  TNtuple* clone_nt = (TNtuple*) nt->CloneTree ();
  clone_nt->SetDirectory (outRootFile);
  clone_nt->Write ();
  outRootFile->Close ();

  inFile->Close ();
  delete inFile;
}
