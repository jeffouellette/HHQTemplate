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

const int NGridX      = 200;   // # of x and y grid cells from hydro output
const int NGridY      = 200; 
const double min_temp = 0.140;
const double time_step_length_fm = 0.0071425; // timestep = (SNAPUPDATE * EPS * AT) / (5.0677 GeV^-1 / fm), e.g. here it is for 141 bins, xmax=5fm

// ==================================================================================
// CREATE HISTOGRAMS TO STORE THE KEY INFORMATION FROM EACH TIMESTEP
const int N_Histograms = 1000;
TH2D* htemperature[N_Histograms+1];
TH2D* henergydensity[N_Histograms+1];
TH2D* hbetax[N_Histograms+1];
TH2D* hbetay[N_Histograms+1];
TH1D* hheavypt[N_Histograms+1];

//Scale factor: the final histogram will have NSCALEFACTOR times more bins in X and Y than the input
const int SCALEFACTOR = 4;

//Current canvas
TCanvas *c_c;

//Histogram with coarse binning
TH2D *h_c;

//Histogram with fine binning
TH2D *h_f;
TH2D *h_f_stage1;
TH2D *h_f_stage2;

//Momentum vectors
TArrow *t_vec[141][141];

//Representation of a point in R^2
struct point {
 float x;
 float y;
 float z;
};

//--------------------------------
// Functions
//--------------------------------

float getLinearInterpolation(float x, float x1, float z1, float x2, float z2)
{
 return z1 + (x - x1) * ((z2 - z1) / (x2 - x1));
}


float getBilinearInterpolation(float x, float y, float x1, float y1, float x2, float y2, float z11, float z12, float z21, float z22)
{
 float aux1 = ((x2 - x) / (x2 - x1)) * z11 + ((x - x1) / (x2 - x1)) * z21;
 float aux2 = ((x2 - x) / (x2 - x1)) * z12 + ((x - x1) / (x2 - x1)) * z22;

 return ((y2 - y) / (y2 - y1)) * aux1 + ((y - y1) / (y2 - y1)) * aux2;
}


int getBoundingRowValues(TH2D *&h, int i, int j, point &pHi, point &pLo)
{
 //Find rightmost point
 point xR;
 xR.x = -999;
 xR.y = -999;
 for (int ibinx = i + 1; ibinx <= h->GetNbinsX(); ibinx++)
 {
  if (h->GetBinContent(ibinx, j) != 0)
  {
   xR.x = ibinx;
   xR.y = j;
   xR.z = h->GetBinContent(ibinx, j);
   break;
  }
 }

 //Find leftmost point
 point xL;
 xL.x = -999;
 xL.y = -999;
 for (int ibinx = i - 1; ibinx > 0; ibinx--)
 {
  if (h->GetBinContent(ibinx, j) != 0)
  {
   xL.x = ibinx;
   xL.y = j;
   xL.z = h->GetBinContent(ibinx, j);
   break;
  }
 }

 if (xL.x != -999 && xR.x != -999)
 {
  pHi = xR;
  pLo = xL;
  return 1;
 }

 return 0;
}


int getBoundingColumnValues(TH2D *&h, int i, int j, point &pHi, point &pLo)
{
 //Find topmost point
 point xR;
 xR.x = -999;
 xR.y = -999;
 for (int ibiny = j + 1; ibiny <= h->GetNbinsY(); ibiny++)
 {
  if (h->GetBinContent(i, ibiny) != 0)
  {
   xR.x = i;
   xR.y = ibiny;
   xR.z = h->GetBinContent(i, ibiny);
   break;
  }
 }

 //Find leftmost point
 point xL;
 xL.x = -999;
 xL.y = -999;
 for (int ibiny = j - 1; ibiny > 0; ibiny--)
 {
  if (h->GetBinContent(i, ibiny) != 0)
  {
   xL.x = i;
   xL.y = ibiny;
   xL.z = h->GetBinContent(i, ibiny);
   break;
  }
 }

 if (xL.x != -999 && xR.x != -999)
 {
  pHi = xR;
  pLo = xL;
  return 1;
 }

 return 0;
}


void findClosestValuesBilinear(TH2D *&hOriginal, TH2D *hNew, int i, int j)
{
 point p1x, p2x;
 point p1y, p2y;

 if (getBoundingRowValues(hOriginal, i, j, p1x, p2x) && getBoundingColumnValues(hOriginal, i, j, p1y, p2y))
 {
  float interp = getBilinearInterpolation(hOriginal->GetXaxis()->GetBinCenter(i), hOriginal->GetYaxis()->GetBinCenter(j), p2x.x, p2x.y, p1y.x, p1y.y, p2x.z,     p1y.z, p1x.z, p2y.z);
  hNew->SetBinContent(i, j, interp);
 }
}


void findClosestValues(TH2D *&hOriginal, TH2D *hNew, int i, int j)
{
 //First check if the (i,j) bin at hand is bounded by points with values in the same row or column.
 //In that case, use simple linear interpolation between the two points
 point p1x, p2x;
 point p1y, p2y;

 if (getBoundingRowValues(hOriginal, i, j, p1x, p2x))
 {
  hNew->SetBinContent(i, j, getLinearInterpolation(hOriginal->GetXaxis()->GetBinCenter(i), hOriginal->GetXaxis()->GetBinCenter(p2x.x), p2x.z, hOriginal->GetXaxis()->GetBinCenter(p1x.x), p1x.z));
 }
 
 if (getBoundingColumnValues(hOriginal, i, j, p1y, p2y))
 {
  hNew->SetBinContent(i, j, getLinearInterpolation(hOriginal->GetYaxis()->GetBinCenter(j), hOriginal->GetYaxis()->GetBinCenter(p2y.y), p2y.z, hOriginal->GetYaxis()->GetBinCenter(p1y.y), p1y.z));
 }
}


void smoothenHistogram()
{
 //Get the number of bins along X and Y
 int nBinsX = h_c->GetNbinsX();
 int nBinsY = h_c->GetNbinsY();

 //Create a new histogram with 'SCALEFACTOR'-times as many bins in each direction
 int newNbinsX = nBinsX * SCALEFACTOR;
 int newNbinsY = nBinsY * SCALEFACTOR;

 cout << "Creating new histogram with " << newNbinsX << " x " << newNbinsY << endl;

 int rangeLowX = h_c->GetXaxis()->GetBinCenter(1) - (h_c->GetXaxis()->GetBinWidth(1) / 2.0);
 int rangeHighX = h_c->GetXaxis()->GetBinCenter(nBinsX) + (h_c->GetXaxis()->GetBinWidth(1) / 2.0);
 int rangeLowY = h_c->GetYaxis()->GetBinCenter(1) - (h_c->GetYaxis()->GetBinWidth(1) / 2.0);
 int rangeHighY = h_c->GetYaxis()->GetBinCenter(nBinsY) + (h_c->GetYaxis()->GetBinWidth(1) / 2.0);

 h_f = new TH2D("h_f", ";x;y", newNbinsX, rangeLowX, rangeHighX, newNbinsY, rangeLowY, rangeHighY);

 //Fill this finer histogram with the values from the coarser histogram
 for (int ibinx = 1; ibinx <= h_c->GetNbinsX(); ibinx++)
 {
  for (int ibiny = 1; ibiny <= h_c->GetNbinsY(); ibiny++)
  {
   float xCenter = h_c->GetXaxis()->GetBinCenter(ibinx);
   float yCenter = h_c->GetYaxis()->GetBinCenter(ibiny);
   float zValue  = h_c->GetBinContent(ibinx, ibiny);

   h_f->Fill(xCenter, yCenter, zValue - h_c->GetBinContent(1, 1));
  }
 }

 //Make a copy of the histogram with fine binning
 h_f_stage1 = (TH2D*) h_f->Clone("h_f_stage1");

 //For a given empty bin, find the four closest non-zero bins
 for (int ibinx = 1; ibinx <= newNbinsX; ibinx++)
 {
  for (int ibiny = 1; ibiny <= newNbinsY; ibiny++)
  {
   if (h_f->GetBinContent(ibinx, ibiny) == 0.0) findClosestValues(h_f, h_f_stage1, ibinx, ibiny);
  }
 }

 //Make a copy of the histogram with fine binning
 h_f_stage2 = (TH2D*) h_f_stage1->Clone("h_f_stage2");

 //For a given empty bin, find the four closest non-zero bins
 for (int ibinx = 1; ibinx <= h_f_stage1->GetNbinsX(); ibinx++)
 {
  for (int ibiny = 1; ibiny <= h_f_stage2->GetNbinsX(); ibiny++)
  {
   if (h_f_stage1->GetBinContent(ibinx, ibiny) == 0.0) findClosestValues(h_f_stage1, h_f_stage2, ibinx, ibiny);
  }
 }
}


void SetStyle (TCanvas* c, TH2D* h) {
  c->SetWindowSize(800, 600);
  c->SetFillColor(kBlue + 3);
  c->SetTickx();
  c->SetTicky();
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.2);
  c->SetBottomMargin(0.13);
  h->SetContour(50);
  h->GetZaxis()->SetRangeUser(0.0, 0.37);
  h->GetXaxis()->SetTitle("x [fm]");
  h->GetYaxis()->SetTitle("y [fm]");
  h->GetXaxis()->SetTitleFont(62);
  h->GetYaxis()->SetTitleFont(62);
  h->GetXaxis()->SetLabelFont(62);
  h->GetYaxis()->SetLabelFont(62);
  //h->GetXaxis()->SetRangeUser(-5, 5);
  //h->GetYaxis()->SetRangeUser(-5, 5);
  h->GetXaxis()->SetAxisColor(kWhite);
  h->GetYaxis()->SetAxisColor(kWhite);
  h->GetXaxis()->SetLabelColor(kWhite);
  h->GetYaxis()->SetLabelColor(kWhite);
  h->GetXaxis()->SetTitleColor(kWhite);
  h->GetYaxis()->SetTitleColor(kWhite);
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
 
  h->GetZaxis()->SetTitleOffset(1.75);
  h->GetZaxis()->SetTitleFont(62);
  h->GetZaxis()->SetLabelFont(62);
  h->GetZaxis()->SetAxisColor(kWhite);
  h->GetZaxis()->SetLabelColor(kWhite);
  h->GetZaxis()->SetTitleColor(kWhite);
  h->GetZaxis()->SetTitle("Temperature [GeV]");
  h->GetZaxis()->SetTitleOffset(1.3);

  h->DrawCopy("COLZ");
  return;
}


void make2DContour()
{
 //Remove white spaces by averaging
 for (int i = 1; i <= h_c->GetNbinsX(); i++)
 {
  for (int j = 1; j <= h_c->GetNbinsY(); j++)
  {
   if (h_c->GetBinContent(i, j) != 0) continue;

   h_c->SetBinContent(i, j, h_c->GetBinContent(i - 1, j));
  }
 }

 for (int i = 1; i <= h_c->GetNbinsX(); i++)
 {
  for (int j = 1; j <= h_c->GetNbinsY(); j++)
  {
   if (h_c->GetBinContent(i, j) != 0) continue;
   
   h_c->SetBinContent(i, j, h_c->GetBinContent(i, j + 1));
  }
 }

 //Attenuate bright background with a modified logistic-type function
 for (int ibinx = 1; ibinx <= h_c->GetNbinsX(); ibinx++)
 {
  for (int ibiny = 1; ibiny <= h_c->GetNbinsY(); ibiny++)
  {
   float cont = h_c->GetBinContent(ibinx, ibiny) / (1 + TMath::Exp(-100 * (h_c->GetBinContent(ibinx, ibiny) - 0.14)));
   h_c->SetBinContent(ibinx, ibiny, cont);
  }
 }
 
 //smoothenHistogram(); // don't always do this - diffusion should not be run on an interpolated grid
 SetStyle (c_c, h_c);
}



void Txt2Hist(int indexspecial=0, // index for collision system. 0-100=He^3-Au, 100-200=d-Au, 200-300=p-Au, 300-400=p-Pb.
              string dir="foo", // input directory (SONIC run directory)
              string dir_out="foo", // output directory
              int timesteps_to_run=0, // how many frames to generate
              double arrowTthresh=0.100, // threshold temperature to plot velocity field
              bool noenergydensity=true,
              int eventNum=0) // for plot labelling
{
  double tmin = 0.140;

  TFile* outFile;
  outFile = new TFile ("nagle-hydro.root", "recreate");

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

  int N_timesteps = 0;

  // now create the files (list of files)
  system(Form("rm %s/input_T.dat",dir.c_str()));
  system(Form("rm %s/input_FO.dat",dir.c_str()));
  system(Form("ls %s/Tcontour_*.dat > input_T.dat",dir.c_str())); // Hydro code outputs files one at a time so they are time ordered
  system(Form("ls %s/FOdata_*.dat > input_FO.dat",dir.c_str()));
  if (!noenergydensity)  system(Form("ls %s/Edata_*.dat > input_ED.dat",dir.c_str()));

  TString files = gSystem->GetFromPipe(Form("ls %s/Tcontour_*.dat | wc -l",dir.c_str()));
  int nfiles = atoi(files.Data());
  cout<<"nfiles: "<<nfiles<<endl;
  N_timesteps = nfiles;
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

  TCanvas* c1[1000];

  // ===============================================================
  // LOOP OVER ALL TIME STEPS FROM THE HYDRO OUTPUT
  if (timesteps_to_run == 0) timesteps_to_run = N_timesteps;
  for (int itime=1; itime<=timesteps_to_run; itime++) {

    double xmax, ymax; // set maximum to extrema of grid
    char fooc1name[1000];
    sprintf(fooc1name,"canvas%d",itime-1);
    sprintf(canvas_name,"c_time_%03d",itime);
    sprintf(hist_temp_name,"h_temp_time_%03d",itime);
    sprintf(hist_ed_name,"h_ed_time_%03d",itime);
    sprintf(hist_betax_name,"h_betax_time_%03d",itime);
    sprintf(hist_betay_name,"h_betay_time_%03d",itime);
    sprintf(hist_heavypt_name,"h_heavypt_time_%03d",itime);

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
    if (!noenergydensity)
      file3 >> ED_name;

    string dummyLine;
    
    data1.open(T_name);
    getline(data1, dummyLine); // skips the first line (comments)
    data1 >> xmax; // the first entry is the time so just do this twice to skip it
    data1 >> xmax;
    data1 >> ymax;
    if (xmax < 0) xmax *= -1;
    if (ymax < 0) ymax *= -1;
    data1.close();
    data1.open(T_name);
    getline(data1, dummyLine);

    data2.open(FO_name);
    
    if (!noenergydensity)
      data3.open(ED_name);
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

    htemperature[itime-1] = new TH2D(hist_temp_name, timeString.Data(), NGridX, -xmax, xmax, NGridY, -ymax, ymax);
    henergydensity[itime-1] = new TH2D(hist_ed_name, timeString.Data(), NGridX, -xmax, xmax, NGridY, -ymax, ymax);
    hbetax[itime-1] = new TH2D(hist_betax_name, timeString.Data(), NGridX, -xmax, xmax, NGridY, -ymax, ymax);
    hbetay[itime-1] = new TH2D(hist_betay_name, timeString.Data(), NGridX, -xmax, xmax, NGridY, -ymax, ymax);
    hheavypt[itime-1] = new TH1D(hist_heavypt_name, timeString.Data(), 25, 0.0, 5.0);
    
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

    // 8/16/2018 J. Ouellette - all formatting now taken care of by Javier's formatting code (see above)

    // set the temperature scale based on the maximum T in the first time step
    c1[itime-1]->cd();
    c1[itime-1]->Draw();

    c_c = c1[itime-1];
    h_c = htemperature[itime-1];
    make2DContour ();

    //char fooname2[1000];
    ////sprintf(fooname2,"Viscous Hydrodynamics, time = %4.3f",1.0+0.0025*100.0*(itime-1));
    //sprintf(fooname2,"Viscous Hydrodynamics, time = %4.3f",time);
    //TText *t1 = new TText(-0.9*xmax,1.02*xmax,fooname2);
    //t1->Draw("same");

    TLatex *tinitial = new TLatex();
    tinitial->SetTextFont (62);
    tinitial->SetTextColor (kWhite);
    if (itime <= 50) {
      TLatex *tinitial = new TLatex();
      tinitial->SetTextFont (62);
      tinitial->SetTextColor (kWhite);
      tinitial->DrawLatexNDC (0.2, 0.8, Form("Event %i", eventNum));
      //      TText *tinitial = new TText(-7.5,8.0,"Cu+Cu Perip. Initial Condition");
      //      TText *tinitial = new TText(-7.5,8.0,"d+Au Central Initial Condition");
      if (indexspecial < 100) {
        tinitial->DrawLatexNDC (0.2, 0.86, "He^{3}+Au Initial Condition");
      } else if (indexspecial < 200) {
        tinitial->DrawLatexNDC (0.2, 0.86, "d+Au Initial Condition");
      } else if (indexspecial < 300) {
        tinitial->DrawLatexNDC (0.2, 0.86, "p+Au Initial Condition");
      } else if (indexspecial < 400) {
        tinitial->DrawLatexNDC (0.2, 0.86, "p+Pb Initial Condition");
      }
    }
    tinitial->DrawLatexNDC (0.2, 0.2, Form("t = %1.2f fm", time));
    delete tinitial;

    data1.close();

    // READ IN THE FREEZE-OUT (FO) TO GET THE FLUID VELOCITIES
    for(int ybin=1; ybin<=NGridY; ybin++) {
      for(int xbin=1; xbin<=NGridX; xbin++) {     
        double beta_x, beta_y;
        data2>>temperature;
        data2>>boost_x;//boost velosity Ux
        data2>>boost_y;//boost velosity Uy
        data2>>extra;
        data2>>extra;
        data2>>extra;
        beta_x = boost_x/ sqrt(1.0+pow(boost_x,2)+pow(boost_y,2));
        beta_y = boost_y/ sqrt(1.0+pow(boost_x,2)+pow(boost_y,2));

        hbetax[itime-1]->SetBinContent(xbin,ybin,beta_x);
        hbetay[itime-1]->SetBinContent(xbin,ybin,beta_y);

        // DRAW BOOST VECTOR FOR EVERY 10TH BIN AND SAVE INTO THE CANVAS !!!
        if(xbin%20==0 && ybin%20==0) {
          double xbincenter = htemperature[itime-1]->GetXaxis()->GetBinCenter(xbin);
          double ybincenter = htemperature[itime-1]->GetYaxis()->GetBinCenter(ybin);
          TArrow *boost_vector = new TArrow (xbincenter,ybincenter,xbincenter+beta_x,ybincenter+beta_y,0.06,"|>");
          boost_vector->SetLineWidth (2);
          boost_vector->SetLineColor (10);
          boost_vector->SetFillColor (10);
          boost_vector->SetFillStyle (1001);
          c1[itime-1]->cd();
          // how about only draw is the cell has a temperature above a certain value...
          // also do not draw arrow head if the vector has length zero (looks bad)
          if ( htemperature[itime-1]->GetBinContent(xbin,ybin) > arrowTthresh &&
               sqrt(beta_x*beta_x + beta_y*beta_y) > 0.12 )          
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

    data2.close();
    if (!noenergydensity)
      data3.close();

    outFile->cd();
    
    htemperature[itime-1]->Write();
    if (!noenergydensity)
      henergydensity[itime-1]->Write();
    hbetax[itime-1]->Write();
    hbetay[itime-1]->Write();
    hheavypt[itime-1]->Write();
    c1[itime-1]->Write();

    if (h_f) delete h_f;
    if (h_f_stage1) delete h_f_stage1;
    if (h_f_stage2) delete h_f_stage2;
    
    
  } // end loop over timesteps

  file1.close();
  file2.close();
  if (!noenergydensity)
    file3.close();

  // done with files
  cout << "Done with reading in HYDRODYNAMIC INFORMATION" << endl;

  outFile->Close();

}
