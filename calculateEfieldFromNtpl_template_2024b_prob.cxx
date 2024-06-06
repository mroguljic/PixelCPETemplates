//! PGM to read Urs' pixel tree and produce summary files.
//! Update to version r23 of the pixel tree
//! Make run specific output files
//! Add hit position to angle output info
//! Require barrel layer one on all tracks
//! Write out hit position and angles
//! Loop over several input files
//! Update to version r28 (not including simhit info)
//! Modify to handle multiple vertices
//! Fix vertex selection bug (3-Dec-2011)
//! Keep 3 layers for each subdetector, add layer info to output file
//! Test 2011 end of year templates
//! Test 2012 June templates
//! Do different FB templates
//! Do full ring and layer templates
//! Use 2016 Summer production templates
//! Separate the new modules by layer
//! Add maximum length and maximum x-size
//! Remove bowing corrections from track angle
//! Add upper/lower row [bias dot] halves
//! Do L1 in new and old FB
//! Update to 2022 IOV4 templates
//! Divide L1 into flipped and unflipped modules
//! For 400V operation of L1
//! Separate flipped and unflipped modules in L1 [unfl = B, flip = F]
//! rH = 1.02
//! use most recent template calibration [early 2024, v1]



#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include "boost/multi_array.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <sys/time.h>
#include "SiPixelTemplate.cc"
static int theVerboseLevel = {0};
#ifdef __arm64__
#include "VVIObj.cc"
#else
#include "VVIObjF.cc"
#endif
#include "SiPixelTemplateReco.cc"

using namespace std;

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"

#define PVMAX 100 
#define MUMAX 100 
#define TRACKMAX 10000 
#define CLPERTRACKMAX 20 
#define CLUSTERMAX 20000 
#define DGPERCLMAX 100  
#define TKPERCLMAX 100  
#define DIGIMAX 100000 

// Global definitions 

// Main program  

int main(int argc, char *argv[])
{
   
   // Local variables 
   bool ydouble[TYSIZE], xdouble[TXSIZE];
   static float xrec, yrec, sigmax, sigmay, probx, proby, probQ;  
   float ddx, ddy;
   float locBz, locBx;
   static int runnumber, nrun, sfile, efile;
   int i, j, ierr, qbin, ID, ptmodule;
   double proba, probmin, probXY=0.;
   static int speed, xhalf, flipped, current_run, external;
   static int tladp1[4], qlad[4]={3, 7, 11, 16};
   static int bptmp[4][8], bpnew[4];
   float cluster[TXSIZE][TYSIZE];
   const int mrow = TXSIZE, mcol = TYSIZE;    
   const int txm1 = TXSIZE - 1, tym1 = TYSIZE - 1;
   //	const float xpitch = 100., ypitch = 150.;
   //	int random(void);
   
   FILE *ifp;
   
   struct timeval now0, now1;
   struct timezone timz;
   long deltas, deltaus;
   double deltat;
   
   printf("Enter lower(0: 0-79)/upper(1: 80-159)/both (2) \n");
   scanf("%d", &xhalf);
   printf("xhalf = %d \n", xhalf);
   
   printf("Enter unflipped(0)/flipped(1)/both(2) \n");
   scanf("%d", &flipped);
   printf("flipped = %d \n", flipped);
   
   
   printf("Enter minimum probability \n");
   scanf("%lf", &probmin);
   printf("require probXY > %lf \n", probmin);
   
   /*  Determine current time */
   
   gettimeofday(&now1, &timz);
   deltas = now1.tv_sec - now0.tv_sec;
   deltaus = now1.tv_usec - now0.tv_usec;
   deltat = ((double)deltaus)/1000000.;
   deltat += (double)deltas;
   gStyle->SetOptStat(0);
   //  gStyle->SetOptFit(10);
   gStyle->SetOptFit(0000);
   gStyle->SetPalette(1,0);
   
   //  The results depend upon temperature and Hall Factor
   // Updated for 2017 Phase 1 and temperatures
   
   double temp[8] = {263.,263.,263.,263.,263.,263.,263.,263.};
//   double rH = 1.11;
   double rH[8] = {1.02,1.02,0.98,0.98,1.02,1.02,1.02,1.02};
   printf("Assumed temperatures: %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf \n",
          temp[0],temp[1],temp[2],temp[3],temp[4],temp[5],temp[6],temp[7]);
   printf("Assumed Hall Factors %lf %lf %lf %lf %lf %lf %lf %lf \n", 
   rH[0],rH[1],rH[2],rH[3],rH[4],rH[5],rH[6],rH[7]);
   
   //  Define Run number and range of subfiles
   
   ifp = fopen("analyze_LT2.txt", "r");
   if (ifp==NULL) {
      printf("no analyze_LT2.txt file found/n");
      return 0;
   }
   
   fscanf(ifp,"%d %d %d %d %d", &runnumber, &nrun, &sfile, &efile, &external);
   fclose(ifp);
   printf("first run number %d, last run number %d, start subfile %d, end subfile %d, external %d\n", runnumber, runnumber+nrun-1, sfile, efile, external);
   
   
   // Initialize template store
   
   speed = -2;
   
   
   // Initialize template store
   
   int bpl1 = 1109;
   int bpl2n = 1098;
   int bpl2 = 1110;
   int bpl3m = 1102;
   int bpl3p = 1103;
   int bpl4m = 1104; 
   int bpl4p = 1105;
   
   bptmp[0][0] = bpl1;
   bptmp[0][1] = bpl1;
   bptmp[0][2] = bpl1;
   bptmp[0][3] = bpl1;
   bptmp[0][4] = bpl1;
   bptmp[0][5] = bpl1;
   bptmp[0][6] = bpl1;
   bptmp[0][7] = bpl1;
   bptmp[1][0] = bpl2;
   bptmp[1][1] = bpl2;
   bptmp[1][2] = bpl2;
   bptmp[1][3] = bpl2;
   bptmp[1][4] = bpl2;
   bptmp[1][5] = bpl2;
   bptmp[1][6] = bpl2;
   bptmp[1][7] = bpl2;
   bptmp[2][0] = bpl3m;
   bptmp[2][1] = bpl3m;
   bptmp[2][2] = bpl3m;
   bptmp[2][3] = bpl3m;
   bptmp[2][4] = bpl3p;
   bptmp[2][5] = bpl3p;
   bptmp[2][6] = bpl3p;
   bptmp[2][7] = bpl3p;
   bptmp[3][0] = bpl4m;
   bptmp[3][1] = bpl4m;
   bptmp[3][2] = bpl4m;
   bptmp[3][3] = bpl4m;
   bptmp[3][4] = bpl4p;
   bptmp[3][5] = bpl4p;
   bptmp[3][6] = bpl4p;
   bptmp[3][7] = bpl4p;
   
   
   for(int lay=0; lay<4; ++lay) {tladp1[lay] = 5*qlad[lay]+1;}
   
   bool newmodule[4][64][8];
   for(int lay=0; lay<4; ++lay) {
      for(int lad=0; lad<64; ++lad) {
         for(int mod=0; mod<8; ++mod) {
            newmodule[lay][lad][mod] = false;
         }
      }
   }
   
   newmodule[1][23][3] = true;
   newmodule[1][21][3] = true;
   newmodule[1][11][3] = true;
   newmodule[1][11][2] = true;
   newmodule[1][11][1] = true;
   newmodule[1][7][4] = true;
   newmodule[1][7][5] = true;
   newmodule[1][7][6] = true;   
   
   bpnew[0] = bpl2n;
   bpnew[1] = bpl2n;
   bpnew[2] = bpl2n;
   bpnew[3] = bpl2n;
   
   // Initialize template store 
   
   std::vector< SiPixelTemplateStore > thePixelTemp_;
   SiPixelTemplate templ(thePixelTemp_);
   
   // Initialize template store 
   
   std::string templates_dir = "templates_dir";
   ID = bpl2;
   templ.pushfile(ID,thePixelTemp_,templates_dir);
   ID = bpl3m;
   templ.pushfile(ID,thePixelTemp_,templates_dir);
   ID = bpl3p;
   templ.pushfile(ID,thePixelTemp_,templates_dir);
   ID = bpl4p;
   templ.pushfile(ID,thePixelTemp_,templates_dir);
   ID = bpl4m;
   templ.pushfile(ID,thePixelTemp_,templates_dir);
   ID = bpl2n;
   templ.pushfile(ID,thePixelTemp_,templates_dir);	
   ID = bpl1;
   templ.pushfile(ID,thePixelTemp_,templates_dir);	   
   ID = bpl1+1;
   templ.pushfile(ID,thePixelTemp_,templates_dir);	       
   
   
   // Selection cuts
   
// Selection cuts

   float pt_cut = 0.5;
   float clusterSizeY_cut[4] = {3.999999,3.999999,3.999999,2.999999};
   float clusterSizeX_cut = 5.0;
//   float residual_cut = 0.005;
   float normChi2_cut = 2.0;
   float clusterCharge_cut = 50.;  //charge in ke per unit thickness
//   float trackQuality_cut = 2.0;
//   int highPurity_cut = 1;
   
   // Histogram and detector parameters
   
   int hist_drift_ = 285;
   int hist_depth_ = 50;
   double min_drift_ = -1000;
   double max_drift_ = 1000;
   double min_depth_ = 000;
   double max_depth_ = 285;
   double thick_ = 0.0285;
   double ypitch_ = 0.0150;
   double hypitch_ = ypitch_/2.;
   double cotbeta_min = clusterSizeY_cut[0]*ypitch_/thick_;
   //   double xpitch_ = 0.0100;
   //   double hxpitch_ = ypitch_/2.;
   
   
   double minfit = 5.;
   double maxfit = 280.;
   //  TF1 *f1 = new TF1("f1","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x",minfit, maxfit);
   TF1 *f1 = new TF1("f1","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x",minfit, maxfit);
   f1->SetParName(0,"offset");
   f1->SetParName(1,"tan#theta_{LA}");
   f1->SetParName(2,"quad term");
   f1->SetParName(3,"cubic term");
   f1->SetParName(4,"quartic term");
   f1->SetParName(5,"quintic term");
   //  f1->SetParName(6,"sextic term");
   //  f1->SetParName(7,"septic term");
   f1->SetParameter(0,0);
   f1->SetParameter(1,0.4);
   f1->SetParameter(2,0.0);
   f1->SetParameter(3,0.0);
   f1->SetParameter(4,0.0);
   f1->SetParameter(5,0.0);
   //  f1->SetParameter(6,0.0);
   //  f1->SetParameter(7,0.0);
   f1->SetLineColor(2);
   f1->SetLineWidth(3);
   
   // Add Linear fit for comparison with previous results
   
   TF1 *f0 = new TF1("f0","[0] + [1]*x",50., 235.);
   f0->SetParName(0,"offset");
   f0->SetParName(1,"tan#theta_{LA}");
   f0->SetParameter(0,0);
   f0->SetParameter(1,0.4);
   f0->SetLineColor(2);
   f0->SetLineWidth(3);
   
   
   ofstream *fAngles[8];
   char fAfile[150];
   sprintf(fAfile,"Efield_output/c_R%6.6d_cotangles_L1B.txt",runnumber);
   fAngles[0] = new ofstream(fAfile, ios::trunc);
   sprintf(fAfile,"Efield_output/c_R%6.6d_cotangles_L1F.txt",runnumber);
   fAngles[1] = new ofstream(fAfile, ios::trunc);
   sprintf(fAfile,"Efield_output/c_R%6.6d_cotangles_L2B.txt",runnumber);
   fAngles[2] = new ofstream(fAfile, ios::trunc);
   sprintf(fAfile,"Efield_output/c_R%6.6d_cotangles_L2N.txt",runnumber);
   fAngles[3] = new ofstream(fAfile, ios::trunc);
   sprintf(fAfile,"Efield_output/c_R%6.6d_cotangles_L3B.txt",runnumber);
   fAngles[4] = new ofstream(fAfile, ios::trunc);
   sprintf(fAfile,"Efield_output/c_R%6.6d_cotangles_L3F.txt",runnumber);
   fAngles[5] = new ofstream(fAfile, ios::trunc);
   sprintf(fAfile,"Efield_output/c_R%6.6d_cotangles_L4B.txt",runnumber);
   fAngles[6] = new ofstream(fAfile, ios::trunc);
   sprintf(fAfile,"Efield_output/c_R%6.6d_cotangles_L4F.txt",runnumber);
   fAngles[7] = new ofstream(fAfile, ios::trunc);
   
	sprintf(fAfile,"Efield_output/lt_bpix%6.6d.txt",runnumber);
	printf("opening file %s\n", fAfile);
	ofstream fClust_bpix(fAfile, ios::trunc );	
   
   
   //  ofstream fAngles("Run2012A_modules1to4_layer3.txt", ios::trunc);
   // 	ofstream fLorentzFit( "lorentzFit.txt", ios::trunc );
   // 	fLorentzFit.precision( 4 );
   // 	fLorentzFit << "module" << "\t" << "layer" << "\t" << "offset" << "\t" << "error" << "\t" << "slope" << "\t" << "error" << "\t" "rel.err" << "\t" "pull" << "\t" << "chi2" << "\t" << "prob" << endl;
   TH2F *h_drift_depth_adc[8], *h_drift_depth_adc2[8], *h_drift_depth_noadc[8];
   TProfile *pq_vs_depth[8];
   
   h_drift_depth_adc[0] = new TH2F("h_drift_depth_adc_1Un", "h_drift_depth_adc_1O",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_adc[1] = new TH2F("h_drift_depth_adc_1Fl", "h_drift_depth_adc_1N",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_adc[2] = new TH2F("h_drift_depth_adc_2B", "h_drift_depth_adc_2B",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_adc[3] = new TH2F("h_drift_depth_adc_2N", "h_drift_depth_adc_2F",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_adc[4] = new TH2F("h_drift_depth_adc_3B", "h_drift_depth_adc_3B",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_adc[5] = new TH2F("h_drift_depth_adc_3F", "h_drift_depth_adc_3F",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_adc[6] = new TH2F("h_drift_depth_adc_4B", "h_drift_depth_adc_4B",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_adc[7] = new TH2F("h_drift_depth_adc_4F", "h_drift_depth_adc_4F",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   
   h_drift_depth_adc2[0] = new TH2F("h_drift_depth_adc2_1Un","h_drift_depth_adc2_1O",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_adc2[1] = new TH2F("h_drift_depth_adc2_1Fl","h_drift_depth_adc2_1N",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_adc2[2] = new TH2F("h_drift_depth_adc2_2B","h_drift_depth_adc2_2B",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_adc2[3] = new TH2F("h_drift_depth_adc2_2N","h_drift_depth_adc2_2F",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_adc2[4] = new TH2F("h_drift_depth_adc2_3B","h_drift_depth_adc2_3B",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_adc2[5] = new TH2F("h_drift_depth_adc2_3F","h_drift_depth_adc2_3F",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_adc2[6] = new TH2F("h_drift_depth_adc2_4B","h_drift_depth_adc2_4B",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_adc2[7] = new TH2F("h_drift_depth_adc2_4F","h_drift_depth_adc2_4F",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   
   h_drift_depth_noadc[0] = new TH2F("h_drift_depth_noadc_1Un",";drift [#mum];depth [#mum]",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_noadc[1] = new TH2F("h_drift_depth_noadc_1Fl",";drift [#mum];depth [#mum]",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_noadc[2] = new TH2F("h_drift_depth_noadc_2B",";drift [#mum];depth [#mum]",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_noadc[3] = new TH2F("h_drift_depth_noadc_2N",";drift [#mum];depth [#mum]",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_noadc[4] = new TH2F("h_drift_depth_noadc_3B",";drift [#mum];depth [#mum]",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_noadc[5] = new TH2F("h_drift_depth_noadc_3F",";drift [#mum];depth [#mum]",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_noadc[6] = new TH2F("h_drift_depth_noadc_4B",";drift [#mum];depth [#mum]",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   h_drift_depth_noadc[7] = new TH2F("h_drift_depth_noadc_4F",";drift [#mum];depth [#mum]",hist_drift_ , min_drift_, max_drift_, hist_depth_, min_depth_, max_depth_);
   
   
   pq_vs_depth[0] = new TProfile("pq_vs_depth_1Un","Avg Pixel Charge_1Un;depth [#mum]",hist_depth_, minfit+10., maxfit-10.,"");
   pq_vs_depth[0]->SetLineColor(2);
   pq_vs_depth[0]->SetMinimum(0.);
   
   pq_vs_depth[1] = new TProfile("pq_vs_depth_1Fl","Avg Pixel Charge_1Fl;depth [#mum]",hist_depth_, minfit+10., maxfit-10.,"");
   pq_vs_depth[1]->SetLineColor(2);
   pq_vs_depth[1]->SetMinimum(0.);
   
   pq_vs_depth[2] = new TProfile("pq_vs_depth_2B","Avg Pixel Charge_2B;depth [#mum]",hist_depth_, minfit+10., maxfit-10.,"");
   pq_vs_depth[2]->SetLineColor(2);
   pq_vs_depth[2]->SetMinimum(0.);
   
   pq_vs_depth[3] = new TProfile("pq_vs_depth_2F","Avg Pixel Charge_2N;depth [#mum]",hist_depth_, minfit+10., maxfit-10.,"");
   pq_vs_depth[3]->SetLineColor(2);
   pq_vs_depth[3]->SetMinimum(0.);
   
   pq_vs_depth[4] = new TProfile("pq_vs_depth_3B","Avg Pixel Charge_3B;depth [#mum]",hist_depth_, minfit+10., maxfit-10.,"");
   pq_vs_depth[4]->SetLineColor(2);
   pq_vs_depth[4]->SetMinimum(0.);
   
   pq_vs_depth[5] = new TProfile("pq_vs_depth_3F","Avg Pixel Charge_3F;depth [#mum]",hist_depth_, minfit+10., maxfit-10.,"");
   pq_vs_depth[5]->SetLineColor(2);
   pq_vs_depth[5]->SetMinimum(0.);
   
   pq_vs_depth[6] = new TProfile("pq_vs_depth_4B","Avg Pixel Charge_4B;depth [#mum]",hist_depth_, minfit+10., maxfit-10.,"");
   pq_vs_depth[6]->SetLineColor(2);
   pq_vs_depth[6]->SetMinimum(0.);
   
   pq_vs_depth[7] = new TProfile("pq_vs_depth_4F","Avg Pixel Charge_4F;depth [#mum]",hist_depth_, minfit+10., maxfit-10.,"");
   pq_vs_depth[7]->SetLineColor(2);
   pq_vs_depth[7]->SetMinimum(0.);
   
   TH1F *h_pt = new TH1F("h_pt","h_pt",100,0,5);
   TH1F *h_ndof = new TH1F("h_ndof","h_ndof",101,0,100);
   TH1I *h_trackQuality = new TH1I("h_trackQuality","h_trackQuality",6,0,6);
   TH1I *h_nHitsPerTrack = new TH1I("h_nHitsPerTrack","h_nHitsPerTrack",50,0,50);
   TH1F *h_qclus = new TH1F("h_qclus","h_qclus",200,0.,50.);
   TH1F * h_chi2 = new TH1F("h_chi2", "h_chi2", 200, 0, 10. );
   TProfile *x_vs_row = new TProfile("x_vs_row","Avg x;row",160, -0.5, 159.,"");
   x_vs_row->SetLineColor(2);
   
   
   
   Int_t run_;
   ULong64_t event_;
   Int_t module_;
   Int_t ladder_;
   Int_t layer_;
   Int_t isflipped_;
   Float_t pt_;
   Float_t p_;
   Float_t eta_;
   Float_t phi_;
   Double_t chi2_;
   Double_t ndof_;
   Int_t trackQuality_;
   Int_t nHitsPerTrack_;
   Int_t isHighPurity_;
   const Int_t maxpix = 500;
   
   struct {
      Int_t npix;
      Float_t row[maxpix];
      Float_t col[maxpix];
      Float_t adc[maxpix];
      Float_t x[maxpix];
      Float_t y[maxpix];
   } pixinfo_;
   
   struct {
      Int_t ncol;
      Int_t dcol[maxpix];
      Float_t adc[maxpix];
      Float_t depth[maxpix];
   } colinfo_;
   
   struct {
      Float_t x;
      Float_t y;
      Double_t alpha;
      Double_t beta;
      Double_t gamma;
   } trackhit_;
   
   struct {
      Float_t x;
      Float_t y;
      Float_t charge;
      Int_t size_x;
      Int_t size_y;
      Int_t maxPixelCol;
      Int_t maxPixelRow;
      Int_t minPixelCol;
      Int_t minPixelRow;
   } clust_;
   
   struct {
      Float_t x;
      Float_t y;
   } rechit_;
   
   current_run = 0;
   
   char infile[150];
   
   int ntclust = 0, nlclust = 0;
   int nclustl1 = 0, nbadcl1 = 0;
   
   for(int runn = runnumber; runn < runnumber+nrun; ++runn) {
      for(int nfile = sfile; nfile <= efile; ++nfile) {
         
         //      if(nfile == 8) continue;
         //      if(nfile == 13) continue;
         //      if(nfile == 25) continue;
         //  Create an input filename for this run
         
         for(int i=0; i<150; ++i) {infile[i] = ' ';}
         if(external == 0) {
            if(nfile < 10) {   
               sprintf(infile,"LorentzTrees/LATree_%6.6d/LA_%6.6d_%1.1d.root",runn,runn,nfile);
            } else {
               if(nfile < 100) { 
                  sprintf(infile,"../LA_calibration_trees/LATree_%6.6d/LA_%6.6d_%2.2d.root",runn,runn,nfile);
               } else {   
                  if(nfile < 1000) { 
                     sprintf(infile,"../LA_calibration_trees/LATree_%6.6d/LA_%6.6d_%3.3d.root",runn,runn,nfile);
                  } else {                  
                     sprintf(infile,"../LA_calibration_trees/LATree_%6.6d/LA_%6.6d_%4.4d.root",runn,runn,nfile);              
                  }               
               }
            }
         } else {
            if(nfile < 10) {   
               sprintf(infile,"../1R23/LA_calibration_trees/LATree_%6.6d/LA_%6.6d_%1.1d.root",runn,runn,nfile);
            } else {
               if(nfile < 100) { 
                  sprintf(infile,"../1R23/LA_calibration_trees/LATree_%6.6d/LA_%6.6d_%2.2d.root",runn,runn,nfile);
               } else {   
                  if(nfile < 1000) { 
                     sprintf(infile,"../1R23/LA_calibration_trees/LATree_%6.6d/LA_%6.6d_%3.3d.root",runn,runn,nfile);
                  } else {                  
                     sprintf(infile,"../1R23/LA_calibration_trees/LATree_%6.6d/LA_%6.6d_%4.4d.root",runn,runn,nfile);              
                  }               
               }
            }         
         }
         printf("file to open is %s\n", infile);
         
         
         TFile *f = new TFile(infile,"READ");
         if(!(f->IsOpen())) {
            printf("file %s not found, skipping\n", infile);
            f->Close();
            continue;
         }
         f->cd();
         
         // fill the histrograms with the ntpl
         TTree * LATree = (TTree*)f->Get("SiPixelLorentzAngleTree_");
         int nentries = LATree->GetEntries();
         LATree->SetBranchAddress("run", &run_);
         LATree->SetBranchAddress("event", &event_);
         LATree->SetBranchAddress("module", &module_);
         LATree->SetBranchAddress("ladder", &ladder_);
         LATree->SetBranchAddress("layer", &layer_);
         LATree->SetBranchAddress("isflipped", &isflipped_);
         LATree->SetBranchAddress("pt", &pt_);
//         LATree->SetBranchAddress("p", &p_);//M
         LATree->SetBranchAddress("eta", &eta_);
         LATree->SetBranchAddress("phi", &phi_);
         LATree->SetBranchAddress("chi2", &chi2_);
         LATree->SetBranchAddress("ndof", &ndof_);
         LATree->SetBranchAddress("trackhit", &trackhit_);
         //      LATree->SetBranchAddress("simhit", &simhit_);
         LATree->SetBranchAddress("npix", &pixinfo_.npix);
         LATree->SetBranchAddress("rowpix", pixinfo_.row);
         LATree->SetBranchAddress("colpix", pixinfo_.col);
         LATree->SetBranchAddress("adc", pixinfo_.adc);
         LATree->SetBranchAddress("xpix", pixinfo_.x);
         LATree->SetBranchAddress("ypix", pixinfo_.y);
         LATree->SetBranchAddress("clust", &clust_); // charge is given in 1000 e
         LATree->SetBranchAddress("rechit", &rechit_);
//         LATree->SetBranchAddress("trackQuality", &trackQuality_);
//         LATree->SetBranchAddress("isHighPurity", &isHighPurity_);
//         LATree->SetBranchAddress("nHitsPerTrack", &nHitsPerTrack_);
         
         cout << "Running over " << nentries << " hits" << endl;
         
         
         for(int ientrie = 0 ; ientrie < nentries; ientrie++){
            //  for(int ientrie = 0 ; ientrie < 20000000; ientrie++){
            LATree->GetEntry(ientrie);
            if(run_ != current_run) {
               current_run = run_;
               printf("current run = %d \n", current_run);
            }
            //         printf("layer/ladder/module = %d/%d/%d \n", layer_, ladder_, module_);
            if(layer_ < 1 || layer_ > 4) continue;
            if(ladder_ < 1 || ladder_ > 64) continue;
            if(module_ < 1 || module_ > 8) continue;

            h_trackQuality->Fill(trackQuality_);
            h_pt->Fill(pt_);
            h_ndof->Fill(ndof_);
            h_nHitsPerTrack->Fill(nHitsPerTrack_);
            double chi2_ndof = chi2_/ndof_;
            h_chi2->Fill(chi2_ndof);
            
            if(pt_< pt_cut)continue;
//            if(chi2_ndof > normChi2_cut)continue;
            //    if(nHitsPerTrack_ < 7) continue;
            //if(ndof_==0)continue;//because of some crashes
            //    if(trackQuality_<trackQuality_cut)continue;
            if(pixinfo_.npix > 25) continue;
            
            bool large_pix = false;
            int hindex = 0;
            for(j=0; j<TXSIZE; ++j) {for(i=0; i<TYSIZE; ++i) {cluster[j][i] = 0.;} }
            for(j=0; j<TXSIZE; ++j) {xdouble[j] = false;}
            for(i=0; i<TYSIZE; ++i) {ydouble[i] = false;}
            //  deltabeta is to remove the bowing corrections to the track angle
//            double deltabeta = 1.07e-3*trackhit_.y;
            double deltabeta = 0.;
            double cotbeta = 1./TMath::Tan(trackhit_.beta-deltabeta);
//            if(fabs(cotbeta) <= cotbeta_min) continue;
            double cotalpha = 1./TMath::Tan(trackhit_.alpha);
            double drdz = sqrt(1.+cotalpha*cotalpha+cotbeta*cotbeta);
            float ccc = clusterCharge_cut*drdz;
            double qclus = clust_.charge/drdz;
            h_qclus->Fill(qclus);
            double drcor = drdz/fabs(cotbeta);
            p_ = pt_*sqrt(1.+cotbeta*cotbeta);
            
            // Find the track ends for later use
            
            colinfo_.ncol = 0;
            int k = 0;
            for (int j = 0; j <  pixinfo_.npix; j++){
               int colpos = static_cast<int>(pixinfo_.col[j]);
               if (pixinfo_.row[j] == 0 || pixinfo_.row[j] == 79 || pixinfo_.row[j] == 80 || pixinfo_.row[j] == 159 || colpos % 52 == 0 || colpos % 52 == 51 ){ large_pix = true;}
            }
            if(large_pix) continue;
            double ravg = (double)(clust_.maxPixelRow + clust_.minPixelRow)/2.;
            x_vs_row->Fill(ravg,(double)rechit_.x);
            
            if(xhalf == 0 && clust_.maxPixelRow > 79) continue;
            if(xhalf == 1 && clust_.minPixelRow < 80) continue;
            
            if(flipped == 0 && isflipped_ == 1) continue;
            if(flipped == 1 && isflipped_ == 0) continue;
            
            
            ID = bptmp[layer_-1][module_-1];
            //         printf("layer/module/ID = %d/%d/%d \n", layer_, module_, ID);
            
            hindex = 2*(layer_-1);
            int lstore = layer_;
//            if(newmodule[layer_-1][ladder_-1][module_-1]) continue;
            if(layer_ == 1 && isflipped_ == 1) { hindex +=1; 
//               ID +=1;
            }
            if(layer_ == 1) {++nclustl1;}
            if(layer_ == 2 && newmodule[layer_-1][ladder_-1][module_-1]){ 
               hindex += 1; 
               ID = bpnew[layer_-1];
               lstore = 5;
            }
            if(layer_ > 2 && module_ > 4) { hindex +=1;}
            if(hindex < 0 || hindex > 8) continue;
            
            //     // is it one of the problematic half ladders? (needs to be excluded)
            //         if( (layer_ == 1 && (ladder_ == 5 || ladder_ == 6 || ladder_ == 15 || ladder_ == 16)) ||(layer_ == 2 && (ladder_ == 8 || ladder_ == 9 || ladder_ == 24 || ladder_ == 25)) ||(layer_ == 3 && (ladder_ ==11 || ladder_ == 12 || ladder_ == 33 || ladder_ == 34)) ) {
            //            continue;
            //         }
       
            float xoff = 10.;
            float yoff = 10.;
            if(pixinfo_.npix > maxpix) continue;           
            

            for (int j = 0; j <  pixinfo_.npix; j++){
               
               // Fill the cluster containers with the cluster
               int ix = pixinfo_.row[j] - clust_.minPixelRow;
               if(ix > txm1) continue;
               int iy = pixinfo_.col[j] - clust_.minPixelCol;
               if(iy > tym1) continue;
               if(pixinfo_.x[j] < xoff) { xoff = pixinfo_.x[j];}
               if(pixinfo_.y[j] < yoff) { yoff = pixinfo_.y[j];}
               cluster[ix][iy] = pixinfo_.adc[j];
               if(pixinfo_.row[j] == 79 || pixinfo_.row[j] == 80) { xdouble[ix] = true;}
               int colpos = static_cast<int>(pixinfo_.col[j]);
               if(colpos % 52 == 0 || colpos % 52 == 51) { ydouble[iy] = true;}
            }
            // Do the template analysis on the cluster 
            SiPixelTemplateReco::ClusMatrix clusterPayload{&cluster[0][0], xdouble, ydouble, mrow,mcol};
            locBx = 1.;
            if(cotbeta < 0.) locBx = -1.;
            locBz = locBx;
            if(cotalpha < 0.) locBz = -locBx;
            ierr = PixelTempReco1D(ID, cotalpha, cotbeta, locBz, locBx, clusterPayload, templ, yrec, sigmay, proby, xrec, sigmax, probx, qbin, speed, probQ);	
            if(ierr != 0) {
               if(ierr != 20) printf("ID %d reco of cotalpha/cotbeta = %f/%f failed with error %d \n", ID, cotalpha, cotbeta, ierr);
               if(layer_ == 1) {++nbadcl1;}
            } else {
               float trackhitx = xrec/10000.+xoff;
               float trackhity = yrec/10000.+yoff;
               ddx = (trackhitx-trackhit_.x)*10000.;
               ddy = (trackhity-trackhit_.y)*10000.;
               proba = probx*proby;
               probXY = proba*(1.-log(proba));
               if(probXY < probmin) continue;
// Don't process/save clusters that are outside the template acceptance
               if((clust_.size_y > 18.) || (clust_.size_x > clusterSizeX_cut)) continue;
               // Save the track angles and momentum; and the module orientation for pixelav simulation
               if(probXY < probmin) continue;                 
               ++ntclust;
               
               *fAngles[hindex] << cotalpha << "\t" << cotbeta << "\t"<< p_ <<"\t"<< isflipped_<<endl;

               if(module_ < 5) {ptmodule = module_ - 5;} else {ptmodule = module_ - 4;}
               fClust_bpix << cotalpha << "  " << cotbeta << "  " << clust_.charge*1000. << "  " << clust_.size_x << "  " << clust_.size_y << "  " << qbin << "  " << probx << "  " << proby << "  "
               << ddx << "  " << ddy << "  " << probQ << "  " << p_ << "  " << isflipped_ << "  " << ptmodule << "  " << lstore << endl;
               
               
               if( (clust_.size_y >= clusterSizeY_cut[layer_-1]) && (clust_.charge < ccc) ){
               
                  ++nlclust;
                  
                  //               if(fabs(trackhity) > 1.6) continue;
                  //               printf("templated hit x/y = %f/%f, rechit x/y = %f/%f \n", trackhitx, trackhity, rechit_.x, rechit_.y);
                  float ylim1 = trackhity - thick_*cotbeta/2.;
                  float ylim2 = trackhity + thick_*cotbeta/2.;
                  float xlim1 = trackhitx - thick_*cotalpha/2.;
                  
                  for (int j = 0; j <  pixinfo_.npix; j++){
                     
                     // Handle the end pixels using the track segment in those pixels
                     
                     float ypixlow = pixinfo_.y[j]-hypitch_;
                     float ypixhigh = pixinfo_.y[j]+hypitch_;
                     if(cotbeta > 0.) {
                        if(ylim1 > ypixlow) ypixlow = ylim1;
                        if(ylim2 < ypixhigh) ypixhigh = ylim2;
                     } else {
                        if(ylim2 > ypixlow) ypixlow = ylim2;
                        if(ylim1 < ypixhigh) ypixhigh = ylim1;
                     }
                     float ypixavg = 0.5*(ypixlow+ypixhigh);
                     float dypix = fabs(ypixhigh-ypixlow);
                     float dx = (pixinfo_.x[j] - xlim1) * 10000.;
                     float dy = (ypixavg - ylim1) * 10000.;
                     float depth = dy * tan(trackhit_.beta);
                     float drift = dx - dy * tan(trackhit_.gamma);
                     double drpix = 10000. * dypix*drcor;
                     //	       double qpix = pixinfo_.adc[j]/drpix;
                     double qpix = pixinfo_.adc[j];
                     h_drift_depth_adc[hindex]->Fill(drift, depth, qpix);
                     h_drift_depth_adc2[hindex]->Fill(drift, depth, qpix*qpix);
                     h_drift_depth_noadc[hindex]->Fill(drift, depth);
                     int dcol = static_cast<int>(pixinfo_.col[j]-clust_.minPixelCol+0.01);
                     if(colinfo_.ncol > maxpix) continue;
                     for(int i = 0; i < colinfo_.ncol; ++i) {
                        if(dcol == colinfo_.dcol[i]) {
                           colinfo_.adc[i] += pixinfo_.adc[j]/drpix;
                           goto skip;
                        }
                     }
                     if(colinfo_.ncol < maxpix) {
                        colinfo_.dcol[colinfo_.ncol] = dcol;
                        colinfo_.depth[colinfo_.ncol] = depth;
                        colinfo_.adc[colinfo_.ncol] = pixinfo_.adc[j]/drpix;
                        ++colinfo_.ncol;
                     }
                  skip: ++k;
                     
                  }
                  for(int i = 0; i < colinfo_.ncol; ++i) {
                     pq_vs_depth[hindex]->Fill(colinfo_.depth[i], colinfo_.adc[i]);
                  }
                  
               }
               
               
            }
         }
         f->Close();
      }
      
   }   
   char sffile[150];
   sprintf(sffile,"Efield_output/c_lorentzsummary%6.6d.txt",runnumber);
   ofstream fLorentzSummary(sffile, ios::trunc );
   
   
   for(int hi = 0; hi < 8; ++hi) {
      TH1F * h_mean = new TH1F("h_mean",";depth [#mum];drift [#mum]", hist_depth_, min_depth_, max_depth_);
      TH1F * h_drift_depth_adc_slice_ = new TH1F("h_slice","h_slice", hist_drift_, min_drift_, max_drift_);
      //loop over bins in depth (z-local-coordinate) (in order to fit slices)
      for( int i = 1; i <= hist_depth_; i++){
         // 				findMean(i, (i_module + (i_layer - 1) * 8));
         double npix = 0;
         
         h_drift_depth_adc_slice_->Reset("ICE");
         
         // determine sigma and sigma^2 of the adc counts and average adc counts
         //loop over bins in drift width
         for( int j = 1; j<= hist_drift_; j++){
            if(h_drift_depth_noadc[hi]->GetBinContent(j, i) >= 1){
               double adc_error2 = (h_drift_depth_adc2[hi]->GetBinContent(j,i) -
                                    h_drift_depth_adc[hi]->GetBinContent(j,i)*h_drift_depth_adc[hi]->GetBinContent(j, i) /
                                    h_drift_depth_noadc[hi]->GetBinContent(j, i)) /  h_drift_depth_noadc[hi]->GetBinContent(j, i);
               h_drift_depth_adc_slice_->SetBinContent(j, h_drift_depth_adc[hi]->GetBinContent(j,i));
               h_drift_depth_adc_slice_->SetBinError(j, sqrt(adc_error2));
               npix += h_drift_depth_noadc[hi]->GetBinContent(j,i);
            }else{
               h_drift_depth_adc_slice_->SetBinContent(j, h_drift_depth_adc[hi]->GetBinContent(j,i));
               h_drift_depth_adc_slice_->SetBinError(j, 0);
            }
         } // end loop over bins in drift width
         
         double mean = h_drift_depth_adc_slice_->GetMean(1);
         double error = 0;
         if(npix != 0){
            error = h_drift_depth_adc_slice_->GetRMS(1) / sqrt(npix);
         }
         
         h_mean->SetBinContent(i, mean);
         h_mean->SetBinError(i, error);
      }// end loop over bins in depth
      
      TH1F *h_mean1 = (TH1F*) h_mean->Clone();
      
      TCanvas * c1 = new TCanvas("c1", "c1", 800, 600);
      
      c1->SaveAs("Efield_output/c1.pdf");
      
      TCanvas * c11 = new TCanvas("c11", "c11", 800, 600);
      
      h_mean->Draw();
      
      h_mean->Fit(f1,"ERQ");
      double p0 = f1->GetParameter(0);
      double e0 = f1->GetParError(0);
      double p1 = f1->GetParameter(1);
      double e1 = f1->GetParError(1);
      double p2 = f1->GetParameter(2);
      double e2 = f1->GetParError(2);
      double p3 = f1->GetParameter(3);
      double e3 = f1->GetParError(3);
      double p4 = f1->GetParameter(4);
      double e4 = f1->GetParError(4);
      double p5 = f1->GetParameter(5);
      double e5 = f1->GetParError(5);
      //   double p6 = f1->GetParameter(6);
      //   double e6 = f1->GetParError(6);
      //   double p7 = f1->GetParameter(7);
      //   double e7 = f1->GetParError(7);
      double chi2 = f1->GetChisquare();
      double prob = f1->GetProb();
      
      char outeps[150];
      char outC[150];
      char label[150];
      bool barrelB;
      int lay;
      if(hi < 8) {
         lay = hi/2;
         if((hi/2*2) == hi) {barrelB = true;} else {barrelB = false;}
      } else {
         lay = hi - 3;
         barrelB = false;
      }
      if(barrelB) {
         sprintf(outeps,"Efield_output/c11_R%6.6d_L%1.1dB.pdf",runnumber,lay+1);
         sprintf(outC,"Efield_output/c11_R%6.6d_L%1.1dB.C",runnumber,lay+1);
         sprintf(label,"barrel_Layer_%1.1dB",lay+1);
      } else {
         sprintf(outeps,"Efield_output/c11_R%6.6d_L%1.1dF.pdf",runnumber,lay+1);
         sprintf(outC,"Efield_output/c11_R%6.6d_L%1.1dF.C",runnumber,lay+1);
         sprintf(label,"barrel_Layer_%1.1dF",lay+1);
      }
      
      c11->SaveAs(outeps);
      c11->SaveAs(outC);
      
      delete c11;
      delete h_mean;
      delete h_drift_depth_adc_slice_;
      cout << "index = " << hi << " " << label << ", quintic fit parameters " << endl;
      
      cout  << p0 << "\t" << e0 << "\t" << p1 << "\t" << e1 << "\t" << e1 / p1 *100. << "\t" << chi2 << "\t" << prob << endl;
      cout  << p2 << "\t" << e2 << "\t" << p3 << "\t" << e3 << "\t" << p4 << "\t" << e4  << "\t" << p5 << "\t" << e5 << endl;
      //   cout  << p6 << "\t" << e6  << "\t" << p7 << "\t" << e7 << endl;
      double cwidth = p1*max_depth_+p2*max_depth_*max_depth_+p3*max_depth_*max_depth_*max_depth_
      +p4*max_depth_*max_depth_*max_depth_*max_depth_
      +p5*max_depth_*max_depth_*max_depth_*max_depth_*max_depth_;
      cout  << "charge width = " << cwidth << endl;
      
      fLorentzSummary << "index = " << hi << " " << label << ", quintic fit parameters " << endl;
      fLorentzSummary << p0 << "\t" << e0 << "\t" << p1 << "\t" << e1 << "\t" << e1 / p1 *100. << "\t" << chi2 << "\t" << prob << endl;
      fLorentzSummary << p2 << "\t" << e2 << "\t" << p3 << "\t" << e3 << "\t" << p4 << "\t" << e4  << "\t" << p5 << "\t" << e5 << endl;
      fLorentzSummary << "charge width = " << cwidth << endl;
      
      double vs = 1.53e9*pow(temp[hi], -0.87);
      double ec = 1.01*pow(temp[hi], 1.55);
      double beta = 2.57e-2*pow(temp[hi], 0.66);
      double ibeta = 1./beta;
      double arg0 = vs*rH[hi]*3.8*1.e-4/ec;
      
      char effile[150];
      if(barrelB) {
         sprintf(effile,"Efield_output/c_lorentzFit%6.6d_L%1.1dB.txt",runnumber,lay+1);
      } else {
         sprintf(effile,"Efield_output/c_lorentzFit%6.6d_L%1.1dF.txt",runnumber,lay+1);
      }
      ofstream fLorentzFit(effile, ios::trunc );
      
      // make electric field profile from the slope
      
      double minplot = minfit+5.;
      double maxplot = maxfit - 4.95;
      double eavg = 0.;
      double neavg = 0.;
      for(double dep = minplot; dep < maxplot ; dep += 5.) {
         double slope = p1+2.*p2*dep+3.*p3*dep*dep+4.*p4*dep*dep*dep+5.*p5*dep*dep*dep*dep;
         //             +6.*p6*dep*dep*dep*dep*dep+7.*p7*dep*dep*dep*dep*dep*dep;
         double arg1 = pow(arg0/slope, beta);
         double efield = ec*pow(arg1-1., ibeta);
         fLorentzFit << dep << "\t" << efield << endl;
         eavg +=efield;
         ++neavg;
      }
      fLorentzFit.close();
      if(barrelB) {
         cout << "Layer " << lay+1 << "B, E_avg = " << eavg/neavg << endl;
      } else {
         cout << "Layer " << lay+1 << "F, E_avg = " << eavg/neavg << endl;
      }
      
      
      TCanvas * c12 = new TCanvas("c12", "c11", 800, 600);
      
      h_mean1->Draw();
      
      h_mean1->Fit(f0,"ERQ");
      p0 = f0->GetParameter(0);
      e0 = f0->GetParError(0);
      p1 = f0->GetParameter(1);
      e1 = f0->GetParError(1);
      
      chi2 = f0->GetChisquare();
      prob = f0->GetProb();
      
      if(barrelB) {
         sprintf(outeps,"Efield_output/c12_R%6.6d_L%1.1dB.pdf",runnumber,lay+1);
      } else {
         sprintf(outeps,"Efield_output/c12_R%6.6d_L%1.1dF.pdf",runnumber,lay+1);
      }
      
      if(barrelB) {
         sprintf(outC,"Efield_output/c12_R%6.6d_L%1.1dB.C",runnumber,lay+1);
      } else {
         sprintf(outC,"Efield_output/c12_R%6.6d_L%1.1dF.C",runnumber,lay+1);
      }
      c12->SaveAs(outeps);
      c12->SaveAs(outC);
      
      delete c12;
      delete h_mean1;
      cout << "index = " << hi << " " << label << ", linear fit parameters " << endl;
      cout  << p0 << "\t" << e0 << "\t" << p1 << "\t" << e1 << "\t" << e1 / p1 *100. << "\t" << chi2 << "\t" << prob << endl << endl;
      fLorentzSummary << "index = " << hi << " " << label << ", linear fit parameters " << endl;
      fLorentzSummary  << p0 << "\t" << e0 << "\t" << p1 << "\t" << e1 << "\t" << e1 / p1 *100. << "\t" << chi2 << "\t" << prob << endl << endl;
      
      
      fAngles[hi]->close();
      
      
   }
   
   printf("total clusters/long clusters = %d/%d \n", ntclust, nlclust);
   printf("total l1 clusters/bad reco L1 = %d/%d \n", nclustl1, nbadcl1);
   
   fLorentzSummary.close();
   
   TCanvas * c2 = new TCanvas("c2", "c2", 800, 600);
   h_pt->Draw("HIST");
   c2->SaveAs("Efield_output/c2.pdf");
   
   TCanvas * c3 = new TCanvas("c3", "c3", 800, 600);
   h_ndof->Draw("HIST");
   c3->SaveAs("Efield_output/c3.pdf");
   
   TCanvas * c4 = new TCanvas("c4", "c4", 800, 600);
   h_nHitsPerTrack->Draw("HIST");
   c4->SaveAs("Efield_output/c4.pdf");
   
   TCanvas * c5 = new TCanvas("c5", "c5", 800, 600);
   h_trackQuality->Draw("HIST");
   c5->SaveAs("Efield_output/c5.pdf");
   
   TCanvas * c7 = new TCanvas("c7", "c7", 800, 600);
   h_qclus->Draw("HIST");
   c7->SaveAs("Efield_output/c7.C");
   c7->SaveAs("Efield_output/c7.pdf");
   
   TCanvas * c8 = new TCanvas("c8", "c8", 800, 600);
   h_chi2->Draw("HIST");
   c8->SaveAs("Efield_output/c8.C");
   c8->SaveAs("Efield_output/c8.pdf");
   
   TCanvas * c13 = new TCanvas("c13", "c13", 800, 600);
   x_vs_row->Draw("HIST");
   c13->SaveAs("Efield_output/c13.pdf");
   
   
   for(int hi = 0; hi < 8; ++hi) {
      
      TCanvas * c6 = new TCanvas("c6", "c6", 800, 600);
      pq_vs_depth[hi]->Draw("HIST");
      char outeps[150];
      bool barrelB;
      int lay;
      //     if(hi < 6) {
      lay = hi/2;
      if((hi/2*2) == hi) {barrelB = true;} else {barrelB = false;}
      //      } else {
      //         lay = hi - 3;
      //         barrelB = false;
      //      }
      if(barrelB) {
         sprintf(outeps,"Efield_output/c6_R%6.6d_L%1.1dB.pdf",runnumber,lay+1);
      } else {
         sprintf(outeps,"Efield_output/c6_R%6.6d_L%1.1dF.pdf",runnumber,lay+1);     
      }
      char outC[150];
      if(barrelB) {
         sprintf(outC,"Efield_output/c6_R%6.6d_L%1.1dB.C",runnumber,lay+1);
      } else {
         sprintf(outC,"Efield_output/c6_R%6.6d_L%1.1dF.C",runnumber,lay+1);
      }
      c6->SaveAs(outeps);
      c6->SaveAs(outC);
      
      //      delete c6;
   }
   
   
   return 0;
} // MAIN__ 

// ***************************************************************** 
//! Calculate 2 gaussianly-distributed random numbers.
// ***************************************************************** 
