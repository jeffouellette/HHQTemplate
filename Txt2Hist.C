#include "TROOT.h"
#include <TSystem.h>
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
#include "TString.h"
#include "TRegexp.h"
#include "TPaletteAxis.h"
#include "TPaveLabel.h"
#include "TPaveText.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// these should be 1 less than the number of grid points given to hydro in params.txt!!
const int NGridX      = 201;   //TODO # of x and y grid cells from hydro output, should be exactly equal to NUMT in params.txt
const int NGridY      = 201;   //TODO
const double min_temp = 0.140;

// ==================================================================================
// CREATE HISTOGRAMS TO STORE THE KEY INFORMATION FROM EACH TIMESTEP
const int N_Histograms = 1000;
TH2D *htemperature[N_Histograms+1];
TH2D *henergydensity[N_Histograms+1];
TH2D *hbetax[N_Histograms+1];
TH2D *hbetay[N_Histograms+1];
TH1D *hheavypt[N_Histograms+1];
TObjArray* outputArray = new TObjArray();



void Txt2Hist(int indexspecial=0, // index for collision system. 0-100=He^3-Au, 100-200=d-Au, 200-300=p-Au, 300-400=p-Pb.
              string dir="foo", // input directory (SONIC run directory)
              string dir_out="foo", // output directory
              int timesteps_to_run=0, // how many frames to generate
              double arrowTthresh=0.100, // threshold temperature to plot velocity field
              bool noenergydensity=true,
              int eventNum=0, // for plot labelling
              double max_temp=0.370) // for scaling the temperature-axis (z-axis)
{
  double tmin = 0.140;

  TFile* outFile = new TFile((dir_out + "/nagle-hydro.root").c_str(), "recreate");
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1);

  // ==================================================================================
  // OPEN THE INPUT HYDRO ASCII FILES 
  // - FIRST THE NAME OF THE FILE THAT LISTS THE TIME STEP FILES
  ifstream file1,file2,file3;
  // - SECOND THE NAME OF THE INDIVIDUAL TIME STEP FILES
  ifstream data1,data2,data3;

  // hydro output files to read in...
  const char *inputfile01="input_T.dat";
  const char *inputfile02="input_FO.dat";
  const char *inputfile03="input_ED.dat"; 

  // now create the files (list of files)
  system(Form("rm %s/input_T.dat",dir.c_str()));
  system(Form("rm %s/input_FO.dat",dir.c_str()));
  system(Form("ls %s/Tcontour_*.dat > input_T.dat",dir.c_str())); // Hydro code outputs files one at a time so they are time ordered
  system(Form("ls %s/FOdata_*.dat > input_FO.dat",dir.c_str()));
  if (!noenergydensity)  system(Form("ls %s/Edata_*.dat > input_ED.dat",dir.c_str()));

  TString files = gSystem->GetFromPipe(Form("ls %s/Tcontour_*.dat | wc -l",dir.c_str()));
  int nfiles = atoi(files.Data());
  cout<<"nfiles: "<<nfiles<<endl;
  int N_timesteps = nfiles;

  file1.open(inputfile01);
  if(!file1) {
    cout<<"No such file found! " << inputfile01 <<endl;
    return;
  }
  file2.open(inputfile02);
  if(!file2) {
    cout<<"No such file found! "<< inputfile02 << endl;
    return;
  }
  if (!noenergydensity) {
    file3.open(inputfile03);
    if(!file3) {
      cout<<"No such file found! "<< inputfile03 << endl;
      return;
    }
  }

  char T_name[200];
  char FO_name[200];
  char ED_name[200];

  char hist_temp_name[200];
  char hist_ed_name[200];
  char hist_betay_name[200];
  char hist_betax_name[200];
  char hist_heavypt_name[200];

  char canvas_name[200];

  double time;
  double temperature;
  double ed;
  double boost_x;
  double boost_y;
  double extra;

  double xpoint, ypoint;

  TCanvas *c1[1000];

  //  TLorentzVector *tquark = new TLorentzVector();

  // ===============================================================
  // LOOP OVER ALL TIME STEPS FROM THE HYDRO OUTPUT
  if (timesteps_to_run == 0) timesteps_to_run = N_timesteps;
  for (int itime=1; itime<=timesteps_to_run; itime++) {

    char fooc1name[1000];
    sprintf(fooc1name,"canvas%d",itime-1);
    c1[itime-1] = new TCanvas(fooc1name,fooc1name,611,53,562,500);
    c1[itime-1]->Range(-12.60663,-12.54642,13.83886,12.49337);
    c1[itime-1]->SetFillColor(0);
    c1[itime-1]->SetBorderMode(0);
    c1[itime-1]->SetBorderSize(2);
    c1[itime-1]->SetRightMargin(0.1451613);
    c1[itime-1]->SetFrameBorderMode(0);
    c1[itime-1]->SetFrameBorderMode(0);
    gPad->SetTicks();

    // READ THE NEXT LINE FROM THE FILELIST NAME AND THEN OPEN THAT SET OF TIME STEP FILES
    file1 >> T_name;
    file2 >> FO_name;
    if (!noenergydensity) file3 >> ED_name;

    string dummyLine;
    
    data1.open(T_name);
    getline(data1, dummyLine); // skips the first line (comments)
    double xmax, ymax; // set maximum to extrema of grid
    data1 >> xmax; // the first entry is the time so just do this twice to skip it
    data1 >> xmax;
    data1 >> ymax;
    if (xmax < 0) xmax *= -1;
    if (ymax < 0) ymax *= -1;
    data1.close();
    data1.open(T_name);
    getline(data1, dummyLine);

    data2.open(FO_name);
    
    if (!noenergydensity)     data3.open(ED_name);
    if(data1 && data2) {
      cout<<"Open files: "<<T_name<<" and "<<FO_name<<" and "<<ED_name<<endl;
    } else if(!data1) {
      cout << "No T_data file found! ==> " << T_name << endl;
      return;
    } else if(!data2) {
      cout << "No FO_data file found! ==> " << FO_name << endl;
      return;
    } 
    if (!noenergydensity) {
      if(!data3) {
	cout << "No ED_data file found! ==> " << ED_name << endl;
	return;
      }
    }
    // Get the time
    TString fileName(T_name);
    TRegexp pattern("[0-9]+\\.[0-9][0-9]+");
    TString timeString("xx.xxx");
    if (fileName.Contains(pattern))
      timeString = fileName(pattern);
    
    //plot the 2D histgram of the temperature with boost velosity
    //plot the 2D histgram of the energy density

    sprintf(canvas_name,"c_time_%03d",itime);
    sprintf(hist_temp_name,"h_temp_time_%03d",itime);
    sprintf(hist_ed_name,"h_ed_time_%03d",itime);
    sprintf(hist_betax_name,"h_betax_time_%03d",itime);
    sprintf(hist_betay_name,"h_betay_time_%03d",itime);
    sprintf(hist_heavypt_name,"h_heavypt_time_%03d",itime);

    htemperature[itime-1] = new TH2D(hist_temp_name,
				     timeString.Data(),
				     NGridX,-xmax,xmax,NGridY,-ymax,ymax);
    henergydensity[itime-1] = new TH2D(hist_ed_name,
				       timeString.Data(),
				       NGridX,-xmax,xmax,NGridY,-ymax,ymax);
    hbetax[itime-1] = new TH2D(hist_betax_name,
			       timeString.Data(),
			       NGridX,-xmax,xmax,NGridY,-ymax,ymax);
    hbetay[itime-1] = new TH2D(hist_betay_name,
			       timeString.Data(),
			       NGridX,-xmax,xmax,NGridY,-ymax,ymax);
    hheavypt[itime-1] = new TH1D(hist_heavypt_name,
				 timeString.Data(),
				 25,0.0,5.0);
    
    // READ IN THE ENERGY DENSITY VALUES AND SET INTO 2D HISTOGRAM FOR THAT TIME STEP
    if (!noenergydensity) {
      double dummy1; double dummy2;
      for(int ybin=1; ybin<=NGridY; ybin++) {
	for(int xbin=1; xbin<=NGridX; xbin++) {  	  
	  data3>>dummy1>>dummy2>>ed;
	  henergydensity[itime-1]->SetBinContent(xbin,ybin,ed);
	}
      }
    }

    // READ IN THE TEMPERATURE VALUES AND SET INTO 2D HISTOGRAM FOR THAT TIME STEP
    /*for(int ybin=1; ybin<=NGridY; ybin++) {
      for(int xbin=1; xbin<=NGridX; xbin++) {
  data1>>xpoint; // fm units
  data1>>ypoint; // fm units
	data1>>temperature; // GEV UNITS
	//	cout << "Temperature values = " << temperature << endl;
//  htemperature[itime-1]->Fill(xpoint,ypoint,temperature);
	htemperature[itime-1]->SetBinContent(htemperature[itime-1]->FindBin(xpoint,ypoint),temperature); //set the bin value with "cell" temperature
      }
    }*/
    time = 1.0+0.0025*100.0*(itime-1);
    for(int ybin=1; ybin<=NGridY; ybin++) {
      for(int xbin=1; xbin<=NGridX; xbin++) {
        data1 >> time;
        data1 >> xpoint;
        data1 >> ypoint;
        data1 >> temperature;
        htemperature[itime-1]->SetBinContent(xbin,ybin,temperature);
      }
    }

    // NOW DRAW THINGS INTO THE CANVAS


    htemperature[itime-1]->SetXTitle("x coordinate [fm]");
    htemperature[itime-1]->SetYTitle("y coordinate [fm]");
    htemperature[itime-1]->SetZTitle("Temperature [GeV]");
    htemperature[itime-1]->GetZaxis()->SetTitleOffset(1.4);
   

    // set the temperature scale based on the maximum T in the first time step
    c1[itime-1]->cd();
    c1[itime-1]->Draw();
    double tmax = htemperature[0]->GetMaximum();
    // should we zero out everything to black below that????????

    // ?????
    cout << "Tmin = " << tmin << " Tmax = " << tmax << endl;
    htemperature[itime-1]->GetZaxis()->SetRangeUser(tmin,tmax);

   TPaletteAxis *palette = new TPaletteAxis(1.01*xmax,-ymax/*-1.05305*xmax*/,1.104265*xmax,ymax/*0.994695*xmax*/,htemperature[itime-1]);
   palette->SetLabelColor(1);
   palette->SetLabelFont(62);
   palette->SetLabelOffset(0.005);
   palette->SetLabelSize(0.06);
   palette->SetTitleOffset(1.4);
   palette->SetTitleSize(0.04);
   palette->SetFillColor(100);
   palette->SetFillStyle(1001);
   htemperature[itime-1]->GetListOfFunctions()->Add(palette,"br");
   htemperature[itime-1]->SetAxisRange(min_temp, max_temp, "Z");
   htemperature[itime-1]->DrawCopy("colz");

    char fooname2[1000];
//    sprintf(fooname2,"Viscous Hydrodynamics, time = %4.3f",1.0+0.0025*100.0*(itime-1));
    sprintf(fooname2,"Viscous Hydrodynamics, time = %4.3f",time);
    TText *t1 = new TText(-0.9*xmax,1.02*xmax,fooname2);
    t1->Draw("same");

    if (itime <= 50) {
      TLatex *tinitial = new TLatex();
      tinitial->DrawLatexNDC(0.2, 0.7, Form("Event %i", eventNum));
      //      TText *tinitial = new TText(-7.5,8.0,"Cu+Cu Perip. Initial Condition");
      //      TText *tinitial = new TText(-7.5,8.0,"d+Au Central Initial Condition");
      tinitial->SetTextColor(2);
      if (indexspecial < 100) {
	//TLatex *tinitial = new TLatex(-7.5,8.0,"He^{3}+Au Central Initial Condition");
	//tinitial->Draw("same");
  tinitial->DrawLatexNDC(0.2, 0.8, "He^{3}+Au Initial Condition");
      } else if (indexspecial < 200) {
	//TLatex *tinitial = new TLatex(-7.5,8.0,"d+Au Central Initial Condition");
	//tinitial->Draw("same");
  tinitial->DrawLatexNDC(0.2, 0.8, "d+Au Initial Condition");
      } else if (indexspecial < 300) {
	//TLatex *tinitial = new TLatex(-7.5,8.0,"p+Au Central Initial Condition");
	//tinitial->Draw("same");
  tinitial->DrawLatexNDC(0.2, 0.8, "p+Au Initial Condition");
      } else if (indexspecial < 400) {
  //TLatex *tinitial = new TLatex(-7.5,8.0,"p+Pb Initial Condition");
  //tinitial->Draw("same");
  tinitial->DrawLatexNDC(0.2, 0.8, "p+Pb Initial Condition");
      }
      delete tinitial;
    
    }

    //    TPaveText *l1 = new  TPaveText(12.0,9.8,11.1,5.0);
    //    l1->AddText("Temperature [GeV]");
    //    l1->SetTextAngle(90.);
    //    l1->Draw("same");


    data1.close();

    // NOW RE-OPEN DATA1 TO GO THROUGH THE TEMPERATURE VALUES
    //data1.open(T_name);
    //getline(data1, dummyLine);
    

    // READ IN THE FREEZE-OUT (FO) TO GET THE FLUID VELOCITIES
    for(int ybin=1; ybin<=NGridY; ybin++) {
      for(int xbin=1; xbin<=NGridX; xbin++) {  	  
//  while(data1) {
//  data1>> xpoint;
//  data1>> ypoint;
//	data1>> temperature; // GEV UNITS
  data2>>temperature;
	data2>>boost_x;//boost velosity Ux
	data2>>boost_y;//boost velosity Uy
	data2>>extra;
	data2>>extra;
	data2>>extra;
	double beta_x = boost_x/ sqrt(1.0+pow(boost_x,2)+pow(boost_y,2));
	double beta_y = boost_y/ sqrt(1.0+pow(boost_x,2)+pow(boost_y,2));

  hbetax[itime-1]->SetBinContent(xbin,ybin,beta_x);
//	hbetax[itime-1]->SetBinContent(hbetax[itime-1]->FindBin(xpoint,ypoint),beta_x);
	hbetay[itime-1]->SetBinContent(xbin,ybin,beta_y);
//  hbetay[itime-1]->SetBinContent(hbetay[itime-1]->FindBin(xpoint,ypoint),beta_y);
//  int xbin = htemperature[itime-1]->GetXaxis()->FindBin(xpoint);
//  int ybin = htemperature[itime-1]->GetYaxis()->FindBin(ypoint);

	// DRAW BOOST VECTOR FOR EVERY 10TH BIN AND SAVE INTO THE CANVAS !!!
	if(xbin%25==0 && ybin%25==0) {
//  if(xbin%3==0 && ybin%3==0) {
	  double xbincenter = htemperature[itime-1]->GetXaxis()->GetBinCenter(xbin);
	  double ybincenter = htemperature[itime-1]->GetYaxis()->GetBinCenter(ybin);
	  TArrow *boost_vector = new TArrow(xbincenter,ybincenter,xbincenter+beta_x,ybincenter+beta_y,0.06,">");
	  boost_vector->SetLineWidth(2);
	  c1[itime-1]->cd();
	  // how about only draw is the cell has a temperature above a certain value...
	  // also do not draw arrow head if the vector has length zero (looks bad)
	  if ( (htemperature[itime-1]->GetBinContent(xbin,ybin) > arrowTthresh) &&
	       (beta_x !=0 || beta_y != 0) )
	       boost_vector->Draw("same");
    	} // end draw boost vector for every 10th bin

      }
    } // end loop over x,y, grid
    char fooname[1000];
    sprintf(fooname,"%s/graphout_%03d.gif", dir_out.c_str(), itime-1);
    c1[itime-1]->SaveAs(fooname);
    //sprintf(fooname,"%s/graphout_%03d.C", dir_out, itime-1);
    //c1[itime-1]->SaveAs(fooname);

    /*
    if (itime == 10) {
      TCanvas *c2 = new TCanvas();
      c2->Divide(2,1);
      c2->cd(1);
      hbetax[itime-1]->Draw("colz");
      c2->cd(2);
      hbetay[itime-1]->Draw("colz");
    }
    */

    //data1.close();
    data2.close();
    if (!noenergydensity)     data3.close();
    
    htemperature[itime-1]->Write();
    henergydensity[itime-1]->Write();
    hbetax[itime-1]->Write();
    hbetay[itime-1]->Write();
    hheavypt[itime-1]->Write();
    c1[itime-1]->Write();

    outputArray->Add(htemperature[itime-1]);
    if (!noenergydensity)     outputArray->Add(henergydensity[itime-1]);
    outputArray->Add(hbetax[itime-1]);
    outputArray->Add(hbetay[itime-1]);
    if (!noenergydensity)     outputArray->Add(hheavypt[itime-1]);
    
  } // end loop over timesteps
  
  file1.close();
  file2.close();
    if (!noenergydensity)   file3.close();

  // done with files
  cout << "Done with reading in HYDRODYNAMIC INFORMATION" << endl;

  outFile->cd();
  outputArray->Write();
  outFile->Close();

}
