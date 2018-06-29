// j.nagle 
//
// 12/21/2012 
// re-edited with a clean version and modified to take input function of temperature for D (not just a/T)
// 01/15/2014
// adding tracking of change in delta-phi versus time

#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TBox.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TArrow.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>

using namespace std;

void SetGraphProps(TGraph* g,
		   Int_t linecolor,
		   Int_t markercolor,
		   Int_t markerstyle,
		   Double_t markersize);

void diffusion(
         const string    inFileName = "nagle-hydro.root", // input nagle-hydro.root file
         const string    inFileName2 = "hydro-b2.0fm.root", // input Ncoll distribution file
         const string    inNcollDistName = "hydro-b2.0fm.root", // input Ncoll distribution name
         const string    inHeavyQuarkPtSpectrumFileName = "heavy_quark_pt.root", // input c, b initial pT spectra file name
         const bool      yesbeauty = false, // only beauty, otherwise only charm
	       const int       NQuarks = 200000000, // 100000
	       int             N_timesteps = 5000, // just give it larger than the 508 steps
	       const TString   output_filename = "fout.root",
         const double    etaOverS = 1,

         const int       eventNum = 1, // parameterA is 3/2pi * etaOverS
	       //double    parameterA = 3.0/(2.*pi), //  / TMath::TwoPi(), // A = DT
	       // for reference, 3/2pi is equivalent to eta/s = 1/4pi
	       // this parameter is overloaded for negative values to use the Dlookup
	       // set parameterA to -1, -2, -3 for temperature dependent cases

	       const double    timeStepScale = 1.0,   // If > 1, skip timesteps, coasting in between. 
	       const bool      usePythiaPt = true,   // option of PYTHIA input quark pT distribution
	       const bool      onepanel = false,      // if you just want the hydro+quark animation
	       const TString   save = "none",         // Canvases: "all", "none", or "last"
	       const double    ScaleDdown = 1.0,      // ignore this 
	       const double    NonInteractTime = 0.0  // how long before you allow the quarks to interact
	       ) 
{

  gROOT->Reset();
  if (gSystem->Getenv("TMPDIR"))
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const double pi = TMath::Pi();
  const double ptMax = 15.0;
  const int numPtBins = 150;
  const int numPhiBins = 120;

  const double QuarkMass                   = (yesbeauty)? 4.20 : 1.27;
  const int    NGridX                      = 200;   // # of x and y grid cells
  const int    NGridY                      = 200; 
  const double time_step_length            = 0.006; 
  const double time_step_length_inverseGeV = time_step_length * (1.0/0.1975);
  const double temperature_cutoff          = 0.1;
  //double etaOverS                          = 1;
  //switch ((int)eventNum/100) {
  //  case 10: etaOverS = 0.5; break;
  //  case 11: etaOverS = 1; break;
  //  case 12: etaOverS = 2; break;
  //  case 13: etaOverS = 3; break;
  //  case 14: etaOverS = 4; break;
  //  default: etaOverS = 0; break;
  //}
  const double parameterA                  = etaOverS * 3.0/(2*pi);

  cout << "DIFFUSION CODE RUNNING with parameter A = " << etaOverS << " x 3/2pi" << endl;

  time_step_length            *= timeStepScale;
  time_step_length_inverseGeV *= timeStepScale;

  TF1 *fkt = new TF1("fkt","TMath::Exp(-x*x/(2.*[0]*[0]))",-10.0,10.0);
  fkt->SetParameter(0,1.0);

  //==============================================================================================
  TH1D *Dlookup; // temperature dependence of D parameter (default as in Teaney and Moore is a/T)
  if (parameterA > 0) {
    Dlookup = new TH1D("Dlookup","Dlookup",100,0.0,0.500); // up to 500 MeV/c
    for (int i=1;i<=100;i++) 
      Dlookup->SetBinContent(i, parameterA / Dlookup->GetBinCenter(i));  // just D=a/T
  } else {
    // read in Dlookup from an external file !!!!
    TFile *inFileDlookup = new TFile("DlookupModel.root");
    if (parameterA == -1) Dlookup  = (TH1D*) inFileDlookup->Get("DlookupA");
    if (parameterA == -2) Dlookup  = (TH1D*) inFileDlookup->Get("DlookupB");
    if (parameterA == -3) Dlookup  = (TH1D*) inFileDlookup->Get("DlookupC");
    if (parameterA == -4) Dlookup  = (TH1D*) inFileDlookup->Get("DlookupD");
    if (parameterA == -5) Dlookup  = (TH1D*) inFileDlookup->Get("DlookupPerfect");
  }
  //==============================================================================================

  //==============================================================================================
  // hydro inputs (good to rename with directories)
  //  TFile *fin2 = new TFile("ncoll-AuAu-b6.5fm.root");         // x,y Ncoll distribution
  //  TFile* inFile  = new TFile("hydro-b6.5fm.root", "read");   // hydro steps: T,E/V,beta,pt
 
  TFile *fin2 = new TFile(inFileName2.c_str(), "read");         // x,y Ncoll distribution
  TFile* inFile = new TFile(inFileName.c_str(), "read");

  //TFile* inFile  = new TFile("hydro-b2.0fm.root", "read");   // hydro steps: T,E/V,beta,pt
  //TFile *fin2 = new TFile(inFileName2.c_str(), "read");         // x,y Ncoll distribution
  //==============================================================================================

  //==============================================================================================
  // input for charm and beauty initial pT shapes
  TFile *fin1 = new TFile(inHeavyQuarkPtSpectrumFileName.c_str(), "read");
  TH1F *quarkpt;

  // The following works in root 6:
  //if (yesbeauty) quarkpt = static_cast <TH1D *> (fin1->Get("b_quark_pt"));
  //else quarkpt = static_cast <TH1D *> (fin1->Get("c_quark_pt"));
  // For root 5 use this:
  if (yesbeauty) quarkpt = (TH1F*)(fin1->Get("b_quark_pt"));
  else quarkpt = (TH1F*)(fin1->Get("charmPtSpectrum"));

  // Alternative to using TH1s from heavy_quark_pt.root
  TF1* ptc = new TF1("ptc", "[0]*x*TMath::Power(x*x + pow([1],2), [2])", 0., 20.);
  ptc->SetParameters(3.3, 2.1, -3.9);    // parameters from arXiv:1205.2396v1 (Cao,Qin, Bass)
  //  ptc->SetParameters(1.0,1.85,-3.52); // just to compare to Teaney and Moore
  ptc->SetLineColor(kRed);

  TF1* ptb = new TF1("ptb", "[0]*x*TMath::Power(x*x + pow([1],2), [2])", 0., 20.);
  ptb->SetParameters(900, 7.5, -4.9);
  ptb->SetLineColor(kBlue);
  //==============================================================================================

  //==============================================================================================
  // Hydro histos for each timestep - to be read from inFile
  const int numhists = N_timesteps+1;
  TH2D *htemperature[numhists];
  TH2D *henergydensity[numhists];
  TH2D *hbetax[numhists];
  TH2D *hbetay[numhists];
  TH1D *hheavypt[numhists];
  TH1D *hheavyraa[numhists];

  TH1D *hheavypt_radius0_5 = new TH1D("hheavypt_radius0_5","hheavypt_radius0_5",numPtBins,0.0,ptMax);
  TH1D *hheavypt_radius1_0 = new TH1D("hheavypt_radius1_0","hheavypt_radius1_0",numPtBins,0.0,ptMax);
  TH1D *hheavypt_radius1_7 = new TH1D("hheavypt_radius1_7","hheavypt_radius1_7",numPtBins,0.0,ptMax);
  TH1D *hheavypt_radius2_5 = new TH1D("hheavypt_radius2_5","hheavypt_radius2_5",numPtBins,0.0,ptMax);

  TH1D *hheavyraa_radius0_5 = new TH1D("hheavyraa_radius0_5","hheavyraa_radius0_5",numPtBins,0.0,ptMax);
  TH1D *hheavyraa_radius1_0 = new TH1D("hheavyraa_radius1_0","hheavyraa_radius1_0",numPtBins,0.0,ptMax);
  TH1D *hheavyraa_radius1_7 = new TH1D("hheavyraa_radius1_7","hheavyraa_radius1_7",numPtBins,0.0,ptMax);
  TH1D *hheavyraa_radius2_5 = new TH1D("hheavyraa_radius2_5","hheavyraa_radius2_5",numPtBins,0.0,ptMax);

  TH1D *hheavypt_radius0_5_orig = new TH1D("hheavypt_radius0_5_orig","hheavypt_radius0_5_orig",numPtBins,0.0,ptMax);
  TH1D *hheavypt_radius1_0_orig = new TH1D("hheavypt_radius1_0_orig","hheavypt_radius1_0_orig",numPtBins,0.0,ptMax);
  TH1D *hheavypt_radius1_7_orig = new TH1D("hheavypt_radius1_7_orig","hheavypt_radius1_7_orig",numPtBins,0.0,ptMax);
  TH1D *hheavypt_radius2_5_orig = new TH1D("hheavypt_radius2_5_orig","hheavypt_radius2_5_orig",numPtBins,0.0,ptMax);

  TH1D *hheavyraa_radius0_5_orig = new TH1D("hheavyraa_radius0_5_orig","hheavyraa_radius0_5_orig",numPtBins,0.0,ptMax);
  TH1D *hheavyraa_radius1_0_orig = new TH1D("hheavyraa_radius1_0_orig","hheavyraa_radius1_0_orig",numPtBins,0.0,ptMax);
  TH1D *hheavyraa_radius1_7_orig = new TH1D("hheavyraa_radius1_7_orig","hheavyraa_radius1_7_orig",numPtBins,0.0,ptMax);
  TH1D *hheavyraa_radius2_5_orig = new TH1D("hheavyraa_radius2_5_orig","hheavyraa_radius2_5_orig",numPtBins,0.0,ptMax);

  TH2D *hheavyraa_radius0_5_particle = new TH2D("hheavyraa_radius0_5_particle","hheavyraa_radius0_5_particle",numPtBins,0.0,ptMax,100,0.0,5.0);
  TH2D *hheavyraa_radius1_0_particle = new TH2D("hheavyraa_radius1_0_particle","hheavyraa_radius1_0_particle",numPtBins,0.0,ptMax,100,0.0,5.0);
  TH2D *hheavyraa_radius1_7_particle = new TH2D("hheavyraa_radius1_7_particle","hheavyraa_radius1_7_particle",numPtBins,0.0,ptMax,100,0.0,5.0);
  TH2D *hheavyraa_radius2_5_particle = new TH2D("hheavyraa_radius2_5_particle","hheavyraa_radius2_5_particle",numPtBins,0.0,ptMax,100,0.0,5.0);

  TH2D *hheavyptphiorig0_5 = new TH2D("hheavyptphiorig0_5","hheavyptphiorig0_5",numPtBins,0.0,ptMax,numPhiBins,-pi,pi);
  TH2D *hheavyptphi0_5 = new TH2D("hheavyptphi0_5","hheavyptphi0_5",numPtBins,0.0,ptMax,numPhiBins,-pi,pi);
  TH2D *hheavyptphiorig1_0 = new TH2D("hheavyptphiorig1_0","hheavyptphiorig1_0",numPtBins,0.0,ptMax,numPhiBins,-pi,pi);
  TH2D *hheavyptphi1_0 = new TH2D("hheavyptphi1_0","hheavyptphi1_0",numPtBins,0.0,ptMax,numPhiBins,-pi,pi);
  TH2D *hheavyptphiorig1_7 = new TH2D("hheavyptphiorig1_7","hheavyptphiorig1_7",numPtBins,0.0,ptMax,numPhiBins,-pi,pi);
  TH2D *hheavyptphi1_7 = new TH2D("hheavyptphi1_7","hheavyptphi1_7",numPtBins,0.0,ptMax,numPhiBins,-pi,pi);
  TH2D *hheavyptphiorig2_5 = new TH2D("hheavyptphiorig2_5","hheavyptphiorig2_5",numPtBins,0.0,ptMax,numPhiBins,-pi,pi);
  TH2D *hheavyptphi2_5 = new TH2D("hheavyptphi2_5","hheavyptphi2_5",numPtBins,0.0,ptMax,numPhiBins,-pi,pi);
  TH2D *hheavyptphiorig = new TH2D("hheavyptphiorig","hheavyptphiorig",numPtBins,0.0,ptMax,numPhiBins,-pi,pi);
  TH2D *hheavyptphi = new TH2D("hheavyptphi","hheavyptphi",numPtBins,0.0,ptMax,numPhiBins,-pi,pi);

  TH1D *hdeltaphiorig = new TH1D("hdeltaphiorig","hdeltaphiorig",30,0.0,pi);
  TH1D *hdeltaphiorig1 = new TH1D("hdeltaphiorig1","hdeltaphiorig1",30,0.0,pi);
  TH1D *hdeltaphiorig2 = new TH1D("hdeltaphiorig2","hdeltaphiorig2",30,0.0,pi);
  TH1D *hdeltaphiorig3 = new TH1D("hdeltaphiorig3","hdeltaphiorig3",30,0.0,pi);
  TH1D *hdeltaphi = new TH1D("hdeltaphi","hdeltaphi",30,0.0,pi);
  hdeltaphi->SetXTitle("charm-anticharm #Delta #phi (rad)");
  TH1D *hdeltaphi1 = new TH1D("hdeltaphi1","hdeltaphi1",30,0.0,pi);
  TH1D *hdeltaphi2 = new TH1D("hdeltaphi2","hdeltaphi2",30,0.0,pi);
  TH1D *hdeltaphi3 = new TH1D("hdeltaphi3","hdeltaphi3",30,0.0,pi);

  TLorentzVector *tquark = new TLorentzVector();
  TLorentzVector *tquark1 = new TLorentzVector();
  TLorentzVector *tquark2 = new TLorentzVector();
  TLatex ltx, ltxn;  ltxn.SetNDC();
  TObjArray* boosts[numhists]; // TArrow* arrays
  TRandom3 ran(0);                  // 0 means new seed every run
  
  TCanvas* cc = new TCanvas("cc", "cc", 10,10,1200,400);
  if (onepanel) cc->SetWindowSize(500,500);  
  if (!onepanel) cc->Divide(3,1);

  for (int itime=0; itime<N_timesteps; itime++) {
    hheavypt[itime]  = new TH1D(Form("h_heavypt_time_%03d", itime+1),Form("h_heavypt_time_%03d", itime+1),numPtBins,0.0,ptMax);
    hheavyraa[itime]  = new TH1D(Form("h_heavyraa_time_%03d", itime+1),Form("h_heavyraa_time_%03d", itime+1),numPtBins,0.0,ptMax);
  }
  // Read in hydro histograms: T, E density, beta, and pt.
  for (int itime=0; itime<N_timesteps; itime++) {
    
    cout << Form("Reading hydro inputs for timestep %d/%d\r", itime+1, N_timesteps);
    
    // check if this information exists in hydro file or not
    if (! inFile->Get(Form("h_temp_time_%03d",    itime+1))) {
      N_timesteps = itime-1;
      cout << "Last hydro timestep = " << N_timesteps << endl;
      break;
    }

    htemperature[itime]   = (TH2D*) inFile->Get(Form("h_temp_time_%03d",    itime+1));
    htemperature[itime]->SetXTitle("x coordinate [fm]");
    htemperature[itime]->SetYTitle("y coordinate [fm]");
    henergydensity[itime] = (TH2D*) inFile->Get(Form("h_ed_time_%03d",      itime+1));
    hbetax[itime]         = (TH2D*) inFile->Get(Form("h_betax_time_%03d",   itime+1));
    hbetay[itime]         = (TH2D*) inFile->Get(Form("h_betay_time_%03d",   itime+1));
    
    // For graphics: boost vectors for every 10th x,y cell
    boosts[itime] = new TObjArray();
    for(int ybin=1; ybin<=NGridY; ybin++) {
      for(int xbin=1; xbin<=NGridX; xbin++) {  	  
	      if(xbin%10==0 && ybin%10==0) {
	        double x = htemperature[itime]->GetXaxis()->GetBinCenter(xbin);
	        double y = htemperature[itime]->GetYaxis()->GetBinCenter(ybin);
	        double beta_x = hbetax[itime]->GetBinContent(xbin,ybin);
	        double beta_y = hbetay[itime]->GetBinContent(xbin,ybin);
	        double length = TMath::Sqrt(beta_x*beta_x + beta_y*beta_y);
	        TArrow *boost_vector = 
	          new TArrow(x, y, x+beta_x, y+beta_y, 0.005*length, ">");
	        boosts[itime]->Add(boost_vector);
	      }
      }
    } // end loop over x,y, grid
  } // end loop over timesteps
  cout << endl;

  //--------------------------------------------------------------------------------------
  // NOW STEP QUARKS THROUGH THE TIME EVOLUTION
  //--------------------------------------------------------------------------------------

  TH1D *hheavypt_orig = new TH1D("hheavypt_orig","hheavypt_orig", 50, 0.0, 5.0);

  // for drawing some example quarks moving around - put in double array of TGraphs
  const int npoints = N_timesteps;
  TGraph *trackheavy[100][npoints];  // map the first 100 for examples
  int maxdraw_trackheavy = 20; // this really draws maxdraw/2 pairs
  for (int it1=0;it1<100;it1++) {
    for (int it2=0;it2<N_timesteps;it2++) {
      trackheavy[it1][it2] = new TGraph();
      trackheavy[it1][it2]->SetMarkerStyle(20);
      trackheavy[it1][it2]->SetMarkerSize(0.4);
      trackheavy[it1][it2]->SetLineColor(it1+1);
      trackheavy[it1][it2]->SetMarkerColor(it1+1);
      trackheavy[it1][it2]->SetLineColor(1);
      trackheavy[it1][it2]->SetMarkerColor(it1%2? kGray : kBlack);
    }
  }

  //==============================================================================================
  // j.nagle - this is how it was from the b6.5 file ??? arggg...
  // 4/18/2018 Jeff Ouellette: static_cast doesn't work right in root 5.. for root 6 this works:
  // .... TH2D *quarkxy = static_cast <TH2D *> (fin2->Get("ncoll_rp"));
  //TH2D* quarkxy= static_cast <TH2D*> (fin2->Get(inNcollDistName.c_str()));
  //
  // for root 5 use this instead:
  TH2D* quarkxy = (TH2D*)(fin2->Get(inNcollDistName.c_str()));
  //==============================================================================================

  if (! quarkpt) cout << "Did not find quark pt distribution TH1D" << endl;
  if (! quarkxy) cout << "Did not find quark x,y distribution TH2D" << endl;

  double qxstore = 0.0;
  double qystore = 0.0;
  double qptstore = 0.0;
  double qphistore = 0.0;
  double qphifinalstore = 0.0;
  double qptfinalstore = 0.0;

  cout << "Staring loop over quarks ..." << endl;
  for(int iquark=1; iquark<=NQuarks; iquark++) {
    
    cout << Form("Quark # %d/%d\r", iquark, NQuarks);
    
    // start with quark initial information
    // in odd numbered events, calculate both momentum for charm and anti-charm
    // and then use the anti-charm for the even numbered event.
    double qpt, qphi, qx, qy;
    if (iquark%2!=0) {

      qx = 0.0; 
      qy = 0.0;
      // Get initial x,y coordinate from the Ncoll spatial distribution
      quarkxy->GetRandom2(qx,qy);

      qpt = 0;
      qphi = ran.Uniform(-TMath::Pi(), TMath::Pi());
      if (usePythiaPt) qpt = quarkpt->GetRandom(); // according to PYTHIA
      else if (yesbeauty) qpt = ptb->GetRandom();
      else qpt = ptc->GetRandom();

      tquark1->SetPxPyPzE( cos(qphi)*qpt, sin(qphi)*qpt , 0.0, sqrt(QuarkMass*QuarkMass + qpt*qpt));
      tquark2->SetPxPyPzE( cos(qphi+TMath::Pi())*qpt, sin(qphi+TMath::Pi())*qpt , 0.0, sqrt(QuarkMass*QuarkMass + qpt*qpt));

      // CONSIDER OPTION OF A KT KICK (SO NOT BACK TO BACK)
      // 1.  RANDOMLY SELECT A KT KICK FROM A GAUSSIAN with sigma = 1 GeV (for example)
      // 2.  in the rest frame of the c/cbar the total mass = 2 x charm quark mass
      // 3.  convert momentum kick in this rest frame into beta boost
      // 4.  then boost each c and cbar using TLorenzVectors

      double kt = fkt->GetRandom();
      // *****
      //      kt = 0.0;
      
      double DiQuarkMass = 2.0 * QuarkMass;
      double ktBetaBoost = TMath::Abs( kt / sqrt(pow(kt,2)+pow(DiQuarkMass,2)) );
      
      double BoostPhi = ran.Uniform(-TMath::Pi(), TMath::Pi());
      tquark1->Boost(ktBetaBoost*sin(BoostPhi),ktBetaBoost*cos(BoostPhi),0.0); 
      tquark2->Boost(ktBetaBoost*sin(BoostPhi),ktBetaBoost*cos(BoostPhi),0.0); 
      // then extract px, py again...
      qpt = tquark1->Pt();
      qphi = tquark1->Phi();// double check....
	
      // store anti-charm information to be grabbed in the next event
      qxstore = qx; // spatial coordinates
      qystore = qy;
      qptstore = tquark2->Pt(); // momentum and phi
      qphistore = tquark2->Phi();

      double deltaphi = tquark1->Phi() - tquark2->Phi();
      if (deltaphi < 0) deltaphi *= -1.0;
      if (deltaphi > pi) deltaphi = 2.0*pi - deltaphi;
      hdeltaphiorig->Fill(deltaphi);
      if (tquark1->Pt()>1.0 && tquark2->Pt()>1.0) hdeltaphiorig1->Fill(deltaphi);
      if (tquark1->Pt()>2.0 && tquark2->Pt()>2.0) hdeltaphiorig2->Fill(deltaphi);
      if (tquark1->Pt()>3.0 && tquark2->Pt()>3.0) hdeltaphiorig3->Fill(deltaphi);

    } else {
      // every other event make the same but opposite direction
      qx = qxstore;
      qy = qystore;
      qpt = qptstore;
      qphi = qphistore;
    }

    // all of the below code uses qx, qy, qpt, qphi... (so need to grab from alternating events above)

    // put into TLorentzVector
    tquark->SetPxPyPzE( cos(qphi)*qpt, sin(qphi)*qpt , 0.0, 
			sqrt(QuarkMass*QuarkMass + qpt*qpt));

    // fill initial value for pT distribution into time step 0
    //----- j.nagle - 12/26/2013 - a little awkward here since kt already included in the original --> NEED TO REMOVE
    hheavypt_orig->Fill(qpt);

    hheavyptphiorig->Fill(qpt, qphi);

    double init_pt;
    // now loop over time steps for diffusion
    for (int istep=0; istep < N_timesteps; istep++) {


      if (istep % (int) timeStepScale) continue;

      if (iquark <= maxdraw_trackheavy) {
	      for (int it1=istep;it1<N_timesteps;it1++) 
	        trackheavy[iquark-1][it1]->SetPoint(istep, qx, qy);
      }
      
      double local_temp  = 0;
      double local_betax = 0;
      double local_betay = 0;

      // from qx,qy coordinate, determine what cell the quark is in... 
      int ibinx = htemperature[istep]->GetXaxis()->FindBin(qx);
      int ibiny = htemperature[istep]->GetYaxis()->FindBin(qy);
      
      // make sure we are inside the grid bounaries
      if (ibinx>=1 && ibinx <=NGridX && ibiny >=1 && ibiny<=NGridY) {
	      local_temp  = htemperature[istep]->GetBinContent(ibinx,ibiny);
	      local_betax = hbetax[istep]->GetBinContent(ibinx,ibiny);
	      local_betay = hbetay[istep]->GetBinContent(ibinx,ibiny);
      }

      if (local_temp >= temperature_cutoff) {
	
	      // boost into rest frame of cell
	      tquark->Boost(-local_betax,-local_betay,0.0);
	      
	      // random kicks
	      //      parameterA = 3.0/(2.*pi);
	      //	double px_kick = ran.Gaus()*sqrt(time_step_length_inverseGeV * 2.0 * pow(local_temp,3) / parameterA);
	      //	double py_kick = ran.Gaus()*sqrt(time_step_length_inverseGeV * 2.0 * pow(local_temp,3) / parameterA);
	      //	double pz_kick = ran.Gaus()*sqrt(time_step_length_inverseGeV * 2.0 * pow(local_temp,3) / parameterA);

	      double Dparameter = Dlookup->GetBinContent(Dlookup->FindBin(local_temp));

	      // if we wanted a momentum dependence to D - take arXiv 1005.0769v1 Figure 20 Riek,Rapp
	      // posit a 0.5 drop in coupling strength in going from 0 GeV to 5 GeV (simple linear for now)
	      bool MomentumDependent = true;
	      if (MomentumDependent) {
	        // D parameter should be increasing with larger pT (i.e. weaker coupling)
	        double scaleValue = 1.0 + 2.0*(tquark->P()/10.0); //  at P=0, sV=1; at P=5, sV = 2.0
	        Dparameter = Dparameter * scaleValue;
	      }

	      double px_kick = ran.Gaus()*sqrt(time_step_length_inverseGeV * 2.0 * pow(local_temp,2) / Dparameter);
	      double py_kick = ran.Gaus()*sqrt(time_step_length_inverseGeV * 2.0 * pow(local_temp,2) / Dparameter);
	      double pz_kick = ran.Gaus()*sqrt(time_step_length_inverseGeV * 2.0 * pow(local_temp,2) / Dparameter);

	      // could rotate and have a different kappa in and out of direction of motion...
	      
	      double quarkbetax = tquark->Px() / tquark->E();
	      double quarkbetay = tquark->Py() / tquark->E();
	      double quarkbetaz = tquark->Pz() / tquark->E();

	      // drag calculation
	      //	double drag_loss_x = time_step_length_inverseGeV * pow(local_temp,2) * quarkbetax / parameterA;
	      //	double drag_loss_y = time_step_length_inverseGeV * pow(local_temp,2) * quarkbetay / parameterA;
	      //	double drag_loss_z = time_step_length_inverseGeV * pow(local_temp,2) * quarkbetaz / parameterA;

	      double drag_loss_x = time_step_length_inverseGeV * pow(local_temp,1) * quarkbetax / Dparameter;
	      double drag_loss_y = time_step_length_inverseGeV * pow(local_temp,1) * quarkbetay / Dparameter;
	      double drag_loss_z = time_step_length_inverseGeV * pow(local_temp,1) * quarkbetaz / Dparameter;

	      // put is all back together 
	      double newpx = tquark->Px() - drag_loss_x + px_kick;
	      double newpy = tquark->Py() - drag_loss_y + py_kick;
	      double newpz = tquark->Pz() - drag_loss_z + pz_kick;

	      // can set a non-interacting time
	      if (((double) istep)*time_step_length <  NonInteractTime) {
	        newpx = tquark->Px();  newpy = tquark->Py();  newpz = tquark->Pz();
	      }

	      tquark->SetPxPyPzE( newpx, newpy , newpz, sqrt(pow(QuarkMass,2)+pow(newpx,2) + pow(newpy,2) + pow(newpz,2)) );
	        
	      // boost back to normal frame (out of fluid rest frame)
	      tquark->Boost(+local_betax, +local_betay, 0.0);
	  
      } // even if below temperature cutoff still fill and just linearly propagate
    
      hheavypt[istep]->Fill(tquark->Pt());
      double radius_orig = sqrt(pow(qxstore,2)+pow(qystore,2));
      double thisphi = tquark->Phi();
      while (TMath::Abs(thisphi) > pi) {
        if (thisphi > 0) thisphi = thisphi - 2*pi;
        else if (thisphi < 0) thisphi = thisphi + 2*pi;
      }
      if (istep == N_timesteps-1) {
	      if (radius_orig > 0.0 && radius_orig < 0.5) {
          hheavypt_radius0_5->Fill(tquark->Pt());
          hheavyraa_radius0_5_particle->Fill(tquark->Pt(),tquark->Pt()/init_pt);
          hheavyptphi0_5->Fill(tquark->Pt(), thisphi);
        } else if (radius_orig > 0.5 && radius_orig < 1.0) {
          hheavypt_radius1_0->Fill(tquark->Pt());
          hheavyraa_radius1_0_particle->Fill(tquark->Pt(),tquark->Pt()/init_pt);
          hheavyptphi1_0->Fill(tquark->Pt(), thisphi);
        } else if (radius_orig > 1.0 && radius_orig < 1.7) {
          hheavypt_radius1_7->Fill(tquark->Pt());
          hheavyraa_radius1_7_particle->Fill(tquark->Pt(),tquark->Pt()/init_pt);
          hheavyptphi1_7->Fill(tquark->Pt(), thisphi);
        } else if (radius_orig > 1.7 && radius_orig < 2.5) {
          hheavypt_radius2_5->Fill(tquark->Pt());
          hheavyraa_radius2_5_particle->Fill(tquark->Pt(),tquark->Pt()/init_pt);
          hheavyptphi2_5->Fill(tquark->Pt(), thisphi);
        }
      } else if (istep == 0) {
        init_pt = tquark->Pt();
	      if (radius_orig > 0.0 && radius_orig < 0.5) {
          hheavypt_radius0_5_orig->Fill(tquark->Pt());
          hheavyptphiorig0_5->Fill(tquark->Pt(), thisphi);
        }
	      else if (radius_orig > 0.5 && radius_orig < 1.0) {
          hheavypt_radius1_0_orig->Fill(tquark->Pt());
          hheavyptphiorig1_0->Fill(tquark->Pt(), thisphi);
        }
	      else if (radius_orig > 1.0 && radius_orig < 1.7) {
          hheavypt_radius1_7_orig->Fill(tquark->Pt());
          hheavyptphiorig1_7->Fill(tquark->Pt(), thisphi);
        }
	      else if (radius_orig > 1.7 && radius_orig < 2.5) {
          hheavypt_radius2_5_orig->Fill(tquark->Pt());
          hheavyptphiorig2_5->Fill(tquark->Pt(), thisphi);
        }
      }
      // update position of quark
      if (tquark->P() > 0) {
	      qx = qx + (tquark->Px()/tquark->E())* time_step_length;
	      qy = qy + (tquark->Py()/tquark->E())* time_step_length;
      } 
      
    } // goto next step

    if (iquark%2!=0) {
      qphifinalstore = tquark->Phi();
      qptfinalstore = tquark->Pt();
    } else {
      double deltaphi = qphifinalstore - tquark->Phi();
      if (deltaphi < 0) deltaphi *= -1.0;
      if (deltaphi > pi) deltaphi = 2.0*pi - deltaphi;
      hdeltaphi->Fill(deltaphi);
      if (tquark->Pt() > 1.0 && qptfinalstore > 1.0) hdeltaphi1->Fill(deltaphi);
      if (tquark->Pt() > 2.0 && qptfinalstore > 2.0) hdeltaphi2->Fill(deltaphi);
      if (tquark->Pt() > 3.0 && qptfinalstore > 3.0) hdeltaphi3->Fill(deltaphi);
    }

    double thisphi = tquark->Phi();
    while (TMath::Abs(thisphi) > pi) {
      if (thisphi > 0) thisphi = thisphi - 2*pi;
      else if (thisphi < 0) thisphi = thisphi + 2*pi;
    }
    hheavyptphi->Fill(tquark->Pt(), thisphi);

  } // end loop over quarks

  //=================================================================================


  cout << "\nDrawing..." << endl;
  
  TH1D *hratio = new TH1D("hratio","hratio",50,0.0,5.0);

  // NOW DRAW THINGS INTO THE CANVAS
  for (int itime=0; itime<N_timesteps; itime++) {

    if (itime % (int) timeStepScale) continue;
    
    // Panel 1: Flow field and quark trajectories
    cc->cd(1);

    Printf("step %d", itime);

    htemperature[itime]->Draw("colz");
    boosts[itime]->Draw();

    for(int iquark=0; iquark<maxdraw_trackheavy; iquark++) {
      if (!trackheavy[iquark][itime])
	Error("", "no trackheavy[%d][%d]", iquark, itime);
      else
	trackheavy[iquark][itime]->Draw("p,same");
    }
    ltx.DrawLatex(-9.0, 10.2, "Temperature");
    
    // Panel 2: pt distribution, original and modified
    if (!onepanel) cc->cd(2);
    hheavypt_orig->SetLineWidth(3);
    hheavypt_orig->SetLineColor(1);
    hheavypt_orig->SetXTitle("Transverse Momentum (GeV/c)");
    hheavypt_orig->SetYTitle("dN/dp_{T} (Initial=black,Current=red)");
    
    // determine the y-axis scale for this figure....
    double maxyrange = hheavypt_orig->GetBinContent(hheavypt_orig->GetMaximumBin());
    if (maxyrange < hheavypt[itime]->GetBinContent(hheavypt[itime]->GetMaximumBin())) 
      maxyrange = hheavypt[itime]->GetBinContent(hheavypt[itime]->GetMaximumBin());
    maxyrange = 1.1 * maxyrange;
    hheavypt_orig->SetMaximum(maxyrange);
    
    if (!onepanel) hheavypt_orig->DrawCopy();
    hheavypt[itime]->SetLineWidth(3);
    hheavypt[itime]->SetLineColor(2);
    if (!onepanel) hheavypt[itime]->DrawCopy("same");
    
    TLine *t1 = new TLine(hheavypt_orig->GetMean(),0.0,hheavypt_orig->GetMean(),10000.0);
    if (!onepanel) t1->Draw("same");
    TLine *t2 = new TLine(hheavypt[itime]->GetMean(),0.0,hheavypt[itime]->GetMean(),10000.0);
    t2->SetLineColor(2);
    if (!onepanel) t2->Draw("same");
    
    if (!onepanel) ltxn.DrawLatex(0.2, 0.9, Form("Time (fm/c): %s", htemperature[itime]->GetTitle()));

    // Panel 3: pt ratio plot ("R_AA")
    if (!onepanel) cc->cd(3);
    for (int ibin=1;ibin<=50;ibin++) {
      if (hheavypt_orig->GetBinContent(ibin) > 0) {
	hratio->SetBinContent(ibin, 
			      hheavypt[itime]->GetBinContent(ibin) / 
			      hheavypt_orig->GetBinContent(ibin));
	hheavyraa[itime]->SetBinContent(ibin, 
			      hheavypt[itime]->GetBinContent(ibin) / 
			      hheavypt_orig->GetBinContent(ibin));
      }
    }

    hratio->SetLineWidth(3);
    hratio->SetLineColor(4);
    hratio->SetMinimum(0.0);
    hratio->SetMaximum(2.0);
    hratio->SetXTitle("Transverse Momentum (GeV/c)");
    hratio->SetYTitle("R_{AA} for this time step");
    maxyrange = hratio->GetBinContent(hratio->GetMaximumBin());
    if (maxyrange < 1.5)
      maxyrange = 1.5;
    maxyrange = 1.1 * maxyrange;
    hratio->SetMaximum(2.0); // try just leaving fixed

    if (!onepanel) hratio->DrawCopy();

    cc->Modified();
    cc->Update();
    cc->Draw();

    // options for saving out information
    if (save.Contains("last")) { 
      if (itime==N_timesteps-(int)timeStepScale) {
	cc->Print(Form("VisualOut/a%.3g.gif",parameterA));
	cc->SaveAs(Form("VisualOut/a%.3g.C",parameterA));
      } else {
	if (itime==N_timesteps) {
	  cc->Print(Form("VisualOut/a%.3g.gif",parameterA));
	  cc->SaveAs(Form("VisualOut/a%.3g.C",parameterA));
	}
      }
    } else if (save.Contains("all")) { 
      cc->Print(Form("VisualOut/canvas_%03d.gif",itime));
    }
    
  } // timesteps

  TCanvas *c2 = new TCanvas();
  c2->cd();
  hdeltaphi1->SetLineWidth(4);
  hdeltaphi1->SetLineColor(1);
  hdeltaphi2->SetLineWidth(4);
  hdeltaphi2->SetLineColor(2);
  hdeltaphi3->SetLineWidth(4);
  hdeltaphi3->SetLineColor(3);
  hdeltaphi->SetMinimum(0.0);
  hdeltaphi->SetLineColor(4);
  hdeltaphi->SetFillColor(5);
  hdeltaphi->Scale(1.0/hdeltaphi->Integral());
  hdeltaphi->SetMinimum(0.0);
  hdeltaphi->DrawCopy("l,f");

  hdeltaphi1->Scale(1.0/hdeltaphi1->Integral());
  hdeltaphi1->DrawCopy("l,same");
  hdeltaphi2->Scale(1.0/hdeltaphi2->Integral());
  hdeltaphi2->DrawCopy("l,same");
  hdeltaphi3->Scale(1.0/hdeltaphi3->Integral());
  hdeltaphi3->DrawCopy("l,same");

  TCanvas *c999 = new TCanvas("c999","c999",10,10,700,700);
  c999->cd();				     
  for (int ibin=1;ibin<=50;ibin++) {
    if (hheavypt_orig->GetBinContent(ibin) > 0) {
      double scale = 1.0;
      scale = hheavypt_radius0_5->Integral()/(hheavypt_radius0_5->Integral()+hheavypt_radius1_0->Integral()+hheavypt_radius1_7->Integral()+hheavypt_radius2_5->Integral());
      hheavyraa_radius0_5->SetBinContent(ibin, 
				       hheavypt_radius0_5->GetBinContent(ibin) / 
				       (scale*hheavypt_orig->GetBinContent(ibin)));
      scale = hheavypt_radius1_0->Integral()/(hheavypt_radius0_5->Integral()+hheavypt_radius1_0->Integral()+hheavypt_radius1_7->Integral()+hheavypt_radius2_5->Integral());
      hheavyraa_radius1_0->SetBinContent(ibin, 
				       hheavypt_radius1_0->GetBinContent(ibin) / 
				       (scale*hheavypt_orig->GetBinContent(ibin)));
      scale = hheavypt_radius1_7->Integral()/(hheavypt_radius0_5->Integral()+hheavypt_radius1_0->Integral()+hheavypt_radius1_7->Integral()+hheavypt_radius2_5->Integral());
      hheavyraa_radius1_7->SetBinContent(ibin, 
				       hheavypt_radius1_7->GetBinContent(ibin) / 
				       (scale*hheavypt_orig->GetBinContent(ibin)));
      scale = hheavypt_radius2_5->Integral()/(hheavypt_radius0_5->Integral()+hheavypt_radius1_0->Integral()+hheavypt_radius1_7->Integral()+hheavypt_radius2_5->Integral());
      hheavyraa_radius2_5->SetBinContent(ibin, 
				       hheavypt_radius2_5->GetBinContent(ibin) / 
				       (scale*hheavypt_orig->GetBinContent(ibin)));
    }
  }
  
  hheavyraa_radius0_5->SetLineWidth(3);
  hheavyraa_radius0_5->SetLineColor(2);
  hheavyraa_radius1_0->SetLineWidth(3);
  hheavyraa_radius1_0->SetLineColor(6);
  hheavyraa_radius1_7->SetLineWidth(3);
  hheavyraa_radius1_7->SetLineColor(4);
  hheavyraa_radius2_5->SetLineWidth(3);
  hheavyraa_radius2_5->SetLineColor(1);
  hheavyraa_radius0_5->SetMinimum(0.0);
  hheavyraa_radius0_5->SetMaximum(2.0);
  hheavyraa_radius0_5->SetXTitle("Transverse Momentum (GeV/c)");
  hheavyraa_radius0_5->SetYTitle("R_{AA}");
  double maxyrange = hheavyraa_radius0_5->GetBinContent(hheavyraa_radius0_5->GetMaximumBin());
  if (maxyrange < 1.5)
    maxyrange = 1.5;
  maxyrange = 1.1 * maxyrange;
  hheavyraa_radius0_5->SetMaximum(2.0); // try just leaving fixed
  hheavyraa_radius0_5->DrawCopy();
  hheavyraa_radius1_0->DrawCopy("same");
  hheavyraa_radius1_7->DrawCopy("same");
  hheavyraa_radius2_5->DrawCopy("same");
  
  // output file
  //==============================================================
  TFile* outFile = new TFile(output_filename,"recreate");

  cc->Write("maincanvas");
  c2->Write("subcanvas");

  for (int i=0;i<N_timesteps;i++) {
    char fooout[100];
    sprintf(fooout,"hheavypt%d",i);
    hheavypt[i]->Write(fooout);
    sprintf(fooout,"hheavyraa%d",i);
    hheavyraa[i]->Write(fooout);
  }
  // write out pT distribution in time steps
  // write out raa in time steps
  // write out canvas of final display ?????
  
  hdeltaphiorig->Write();
  hdeltaphiorig1->Write();
  hdeltaphiorig2->Write();
  hdeltaphiorig3->Write();
  hdeltaphi->Write();
  hdeltaphi1->Write();
  hdeltaphi2->Write();
  hdeltaphi3->Write();

  hheavyptphiorig0_5->Write();
  hheavyptphi0_5->Write();
  hheavyptphiorig1_0->Write();
  hheavyptphi1_0->Write();
  hheavyptphiorig1_7->Write();
  hheavyptphi1_7->Write();
  hheavyptphiorig2_5->Write();
  hheavyptphi2_5->Write();
  hheavyptphiorig->Write();
  hheavyptphi->Write();

  hheavypt_radius0_5_orig->Write();
  hheavypt_radius1_0_orig->Write();
  hheavypt_radius1_7_orig->Write();
  hheavypt_radius2_5_orig->Write();
  hheavypt_radius0_5->Write();
  hheavypt_radius1_0->Write();
  hheavypt_radius1_7->Write();
  hheavypt_radius2_5->Write();

  hheavyraa_radius0_5->Write();
  hheavyraa_radius1_0->Write();
  hheavyraa_radius1_7->Write();
  hheavyraa_radius2_5->Write();
  hheavyraa_radius0_5_particle->Write();
  hheavyraa_radius1_0_particle->Write();
  hheavyraa_radius1_7_particle->Write();
  hheavyraa_radius2_5_particle->Write();

  outFile->Close();

}  

void SetGraphProps(TGraph* g,
		   Int_t linecolor,
		   Int_t markercolor,
		   Int_t markerstyle,
		   Double_t markersize) 
{
  g->SetLineColor(linecolor);
  g->SetMarkerColor(markercolor);
  g->SetMarkerStyle(markerstyle);
  g->SetMarkerSize(markersize);
  g->SetLineWidth(2);
}
