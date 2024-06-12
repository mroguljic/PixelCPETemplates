//! \file template_code6d.cpp
//!
//! Template hit reconstruction algorithms, add FPix to templates, add double pixels to templates, 
//! change dp calling seq (Add PSI46v2 response function and allow for zero-suppressed ROC output)
//! Change to Root historgraams for long-term compatibility
//! Add angle vs resolution for templates and "Standard Algorithms"
//! Tune CMSSW simulation for template 1 reconstruction
//! Change standard algorithm to always use edge method for y-reconstruction
//! Add Estar template number 4
//! Do cosmics
//! Add clustering  algorithm
//! Update to two threshold response function
//! Take simulation charge scale from template
//! replace rechit position with simhit position
//! set module to forward L1
//! Put reconstructed position in trackhit_
//! Update to 2017 versions of everything
//! Use L1 VCAL defs


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
#include "SiPixelGenError.cc"
#include "SiPixelTemplate.cc"
static int theVerboseLevel = {2};
#include "SiPixelTemplateReco.cc"
#include "VVIObj.cc"
#include "PixelGeneric2D.cc"

using namespace std;

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TObject.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TTree.h"

// Global definitions 
// Global definitions 

struct Pixinfo
{
	int npix;
	float row[100];
	float col[100];
	float adc[100];
	float x[100];
	float y[100];
};

struct Hit{
	float x;
	float y;
	double alpha;
	double beta;
	double gamma;
}; 

struct Clust {
	float x;
	float y;
	float charge;
	int size_x;
	int size_y;
	int maxPixelCol;
	int maxPixelRow;
	int minPixelCol;
	int minPixelRow;
};

struct Rechit {
	float x;
	float y;
};


struct TTRechit {
	float x;
	float y;
	float sigmax;
	float sigmay;
	float probx;
	float proby;
	int qbin;
};

// Main program  

int main(int argc, char *argv[])
{
	typedef boost::multi_array<float, 2> array_2d;
	typedef boost::multi_array<bool, 2> barray_2d;
	

    // Local variables 
	std::vector<float> pvec(6), wgauss(TYSIZE), vgauss(TYSIZE), xgauss(TYSIZE), ygauss(TYSIZE), zgauss(TYSIZE);
    float pixin[TXSIZE][TYSIZE];
    static bool fpix;
    bool ydouble[TYSIZE], xdouble[TXSIZE];
    int mrow = TXSIZE, mcol = TYSIZE;	
    static float thick, xsize, ysize, noise, zcen, xcmssw, ycmssw, sxcmssw, sycmssw, gain_frac, q100_frac, common_frac, readout_noise;
	static float xhit, yhit, xrec, yrec, sigmax, sigmay, probx, proby, signal, cotalpha, cotbeta, qclust, probQ, qscale;
	double proba, probXY, probmin;
	static float yreccm, yhitcm, xreccm, xhitcm, xcmsswcm, ycmsswcm;
    static int ndata, nfile, neh, nevent, ID, nbad, non_linear, icol, ndcol; 
	static vector<int> nbin(5,0);
    int i, j, k, ierr, qbin, nypix, nxpix, etabin, jmin, jmax, imin, imax, numadd, npix, idcol;
	int colmin, colmax, col, rowmin, rowmax, row;
	double dx, dy, eta, etas, dxc, dyc, log10proby, log10probx, tote, bade, weight, alpha, beta, gamma, adc;
	static int iyd, ixd, speed;	
	static float q100, q101, qmax; 
	static double gain = 3.19;
	static double ped = 16.46;
	static double p0 = 0.01218;
	static double p1 = 0.711;
	static double p2 = 203.;
	static double p3 = 148.;	
	static double vcal = 50.;	
	static double vcaloffst = 670.;	
   static float qseed = 1500.;
	float qin, yfrac, xfrac, locBz, locBx;
    static char infile[80], header[80], c, outfile0[80], outfile1[80], outfile2[80], ntupfile[80];
	int triplg(std::vector<float>&);
//	int random(void);
	float cluster[TXSIZE][TYSIZE], clust[TXSIZE][TYSIZE], rclust[TXSIZE][TYSIZE];
	bool bclust[TXSIZE][TYSIZE];

    std::pair<int, int> pixel, max;

    FILE *ifp;
	
   struct timeval now0, now1;
   struct timezone timz;
   long deltas, deltaus;
   double deltat;
        
	//  Read which data and inputs to use (use c file i/o which is much faster than c++ i/o) 
	
	ifp = fopen("q_dist_2t.txt", "r");
	if (ifp==NULL) {
      printf("no pixel initialization file found/n");
      return 0;
	}
	
	
	fscanf(ifp,"%d %d %f %f %f %f %f %f %f %d", &ndata, &nfile, &noise, &q100, &q101, &q100_frac, &common_frac, &gain_frac, &readout_noise, &non_linear);
	fclose(ifp);
	printf("data file %d, mc file %d noise = %f, threshold0 = %f, threshold1 = %f, rms threshold frac = %f, common_frac = %f, gain fraction = %f, readout noise = %f, nonlinear_resp = %d \n", ndata, nfile, noise, q100, q101, q100_frac, common_frac, gain_frac, readout_noise, non_linear);

    //  Create an input filename for this run 
    
    if(nfile < 10000) {
        
        sprintf(infile,"simulated_clusters/template_events_d%4.4d.out",nfile);
        
    } else {
        
        sprintf(infile,"simulated_clusters/template_events_d%5.5d.out",nfile);
        
    }

//  Open input file and read header info 

	ifp = fopen(infile, "r");
    if (ifp==NULL) {
      printf("no pixel data file found/n");
      return 0;
    }
	
// Read-in a header string first and print it    
    
    for (i=0; (c=getc(ifp)) != '\n'; ++i) {
       if(i < 79) {header[i] = c;}
    }
	if(i > 78) {i=78;}
	header[i+1] ='\0';
    printf("%s\n", header);
	
	double  halfxs=300.;
	int nx=120;	
	gStyle->SetOptFit(101);
	gStyle->SetHistLineWidth(2);
	static vector<TH1F*> hp(24);
	hp[0] = new TH1F("h101","dy_temp (all sig); #Deltay (#mum)",nx,-halfxs,halfxs);
	hp[1] = new TH1F("h102","dy_temp (signal @> 1.5mn); #Deltay (#mum)",nx,-halfxs,halfxs);      
	hp[2] = new TH1F("h103","dy_temp (1.5mn @> signal @> 1.0mn); #Deltay (#mum)",nx,-halfxs,halfxs);      
	hp[3] = new TH1F("h104","dy_temp (1.0mn @> signal @> 0.85mn); #Deltay (#mum)",nx,-halfxs,halfxs);     
	hp[4] = new TH1F("h105","dy_temp (0.85mn @> signal); #Deltay (#mum)",nx,-halfxs,halfxs);      
	hp[5] = new TH1F("h201","log10(proby) (all sig)",nx,-12.,0.);
	hp[6] = new TH1F("h202","log10(proby) (signal @> 1.5mn)",nx,-12.,0.);      
	hp[7] = new TH1F("h203","log10(proby) (1.5mn @> signal @> 1.0mn)",nx,-12.,0.);      
	hp[8] = new TH1F("h204","log10(proby) (1.0mn @> signal @> 0.85mn)",nx,-12.,0.);     
	hp[9] = new TH1F("h205","log10(proby) (0.85mn @> signal)",nx,-12.,0.);      
	hp[10] = new TH1F("h106","dx_temp (all sig); #Deltax (#mum)",nx,-halfxs,halfxs);
	hp[11] = new TH1F("h107","dx_temp (signal @> 1.5mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
	hp[12] = new TH1F("h108","dx_temp (1.5mn @> signal @> 1.0mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
	hp[13] = new TH1F("h109","dx_temp (1.0mn @> signal @> 0.85mn); #Deltax (#mum)",nx,-halfxs,halfxs);      
	hp[14] = new TH1F("h110","dx_temp (0.85mn @> signal); #Deltax (#mum)",nx,-halfxs,halfxs);    
	hp[15] = new TH1F("h206","log10(probx) (all sig)",nx,-12.,0.);
	hp[16] = new TH1F("h207","log10(probx) (signal @> 1.5mn)",nx,-12.,0.);      
	hp[17] = new TH1F("h208","log10(probx) (1.5mn @> signal @> 1.0mn)",nx,-12.,0.);      
	hp[18] = new TH1F("h209","log10(probx) (1.0mn @> signal @> 0.85mn)",nx,-12.,0.);     
	hp[19] = new TH1F("h210","log10(probx) (0.85mn @> signal)",nx,-12.,0.);      
	hp[20] = new TH1F("h300","cotbeta (probx<10-3)",nx,-10.,10.);      
	hp[21] = new TH1F("h301","cotalpha (probx<10-3)",nx,-0.25,0.25);   
	hp[22] = new TH1F("h401","dy_cmssw (all sig); #Deltay (#mum)",nx,-halfxs,halfxs);
	hp[23] = new TH1F("h406","dx_cmssw (all sig); #Deltax (#mum)",nx,-halfxs,halfxs);
	
// Set style for the the histograms	
	
	for(i=0; i<24; ++i) {
	   hp[i]->SetLineColor(2);
	   hp[i]->SetFillColor(38);
	}

// Make some profile Histograms

	static vector<TProfile*> pp(48);
	pp[0] = new TProfile("pqhytmp","1.5>Q/Q_avg>1, y_temp",11,0,2.75,"s"); 
	pp[1] = new TProfile("pqhycms","1.5>Q/Q_avg>1, y_cmssw",11,0,2.75,"s"); 
	pp[2] = new TProfile("pqlytmp","1>Q/Q_avg, y_temp",11,0,2.75,"s"); 
	pp[3] = new TProfile("pqlycms","1>Q/Q_avg, y_cmssw",11,0,2.75,"s"); 
	pp[4] = new TProfile("pqhxtmp","1.5>Q/Q_avg>1, x_temp",11,0,2.75,"s"); 
	pp[5] = new TProfile("pqhxcms","1.5>Q/Q_avg>1, x_cmssw",11,0,2.75,"s"); 
	pp[6] = new TProfile("pqlxtmp","1>Q/Q_avg, x_temp",11,0,2.75,"s"); 
	pp[7] = new TProfile("pqlxcms","1>Q/Q_avg, x_cmssw",11,0,2.75,"s"); 
	pp[8] = new TProfile("pycfrac2","dyc vs yfrac, 2 pix clust; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[9] = new TProfile("pycfrac3","dyc vs yfrac, 3 pix clust; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[10] = new TProfile("pycfrac4","dyc vs yfrac, 4 pix clust; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[11] = new TProfile("pxcfrac2","dxc vs xfrac, 2 pix clust; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[12] = new TProfile("pxcfrac3","dxc vs xfrac, 3 pix clust; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[13] = new TProfile("pxcfrac4","dxc vs xfrac, 4 pix clust; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[14] = new TProfile("pyfrac2","dy vs yfrac, 2 pix clust; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[15] = new TProfile("pyfrac3","dy vs yfrac, 3 pix clust; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[16] = new TProfile("pyfrac4","dy vs yfrac, 4 pix clust; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[17] = new TProfile("pxfrac2","dx vs xfrac, 2 pix clust; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[18] = new TProfile("pxfrac3","dx vs xfrac, 3 pix clust; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[19] = new TProfile("pxfrac4","dx vs xfrac, 4 pix clust; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[20] = new TProfile("pycfrac2q1","dyc vs yfrac, 2 pix clust, 1.5>Q/Q_avg>1; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[21] = new TProfile("pycfrac3q1","dyc vs yfrac, 3 pix clust, 1.5>Q/Q_avg>1; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[22] = new TProfile("pycfrac4q1","dyc vs yfrac, 4 pix clust, 1.5>Q/Q_avg>1; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[23] = new TProfile("pxcfrac2q1","dxc vs xfrac, 2 pix clust, 1.5>Q/Q_avg>1; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[24] = new TProfile("pxcfrac3q1","dxc vs xfrac, 3 pix clust, 1.5>Q/Q_avg>1; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[25] = new TProfile("pxcfrac4q1","dxc vs xfrac, 4 pix clust, 1.5>Q/Q_avg>1; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[26] = new TProfile("pyfrac2q1","dy vs yfrac, 2 pix clust, 1.5>Q/Q_avg>1; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[27] = new TProfile("pyfrac3q1","dy vs yfrac, 3 pix clust, 1.5>Q/Q_avg>1; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[28] = new TProfile("pyfrac4q1","dy vs yfrac, 4 pix clust, 1.5>Q/Q_avg>1; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[29] = new TProfile("pxfrac2q1","dx vs xfrac, 2 pix clust, 1.5>Q/Q_avg>1; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[30] = new TProfile("pxfrac3q1","dx vs xfrac, 3 pix clust, 1.5>Q/Q_avg>1; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[31] = new TProfile("pxfrac4q1","dx vs xfrac, 4 pix clust, 1.5>Q/Q_avg>1; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[32] = new TProfile("pycfrac2q2","dyc vs yfrac, 2 pix clust, 1>Q/Q_avg; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[33] = new TProfile("pycfrac3q2","dyc vs yfrac, 3 pix clust, 1>Q/Q_avg; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[34] = new TProfile("pycfrac4q2","dyc vs yfrac, 4 pix clust, 1>Q/Q_avg; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[35] = new TProfile("pxcfrac2q2","dxc vs xfrac, 2 pix clust, 1>Q/Q_avg; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[36] = new TProfile("pxcfrac3q2","dxc vs xfrac, 3 pix clust, 1>Q/Q_avg; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[37] = new TProfile("pxcfrac4q2","dxc vs xfrac, 4 pix clust, 1>Q/Q_avg; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[38] = new TProfile("pyfrac2q2","dy vs yfrac, 2 pix clust, 1>Q/Q_avg; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[39] = new TProfile("pyfrac3q2","dy vs yfrac, 3 pix clust, 1>Q/Q_avg; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[40] = new TProfile("pyfrac4q2","dy vs yfrac, 4 pix clust, 1>Q/Q_avg; y_{frac}; #Deltay (#mum)",20,0.,1.0," "); 
	pp[41] = new TProfile("pxfrac2q2","dx vs xfrac, 2 pix clust, 1>Q/Q_avg; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[42] = new TProfile("pxfrac3q2","dx vs xfrac, 3 pix clust, 1>Q/Q_avg; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[43] = new TProfile("pxfrac4q2","dx vs xfrac, 4 pix clust, 1>Q/Q_avg; x_{frac}; #Deltax (#mum)",20,0.,1.0," "); 
	pp[44] = new TProfile("pqhxtmpvsa","1.5>Q/Q_avg>1, x_temp",25,-1.125,1.125,"s"); 
	pp[45] = new TProfile("pqhxcmsvsa","1.5>Q/Q_avg>1, x_cmssw",25,-1.125,1.125,"s"); 
	pp[46] = new TProfile("pqlxtmpvsa","1>Q/Q_avg, x_temp",25,-1.125,1.125,"s"); 
	pp[47] = new TProfile("pqlxcmsvsa","1>Q/Q_avg, x_cmssw",25,-1.125,1.125,"s"); 

// Set style for the the profile histograms	
	
	for(i=0; i<6; ++i) {
	   pp[i+8]->SetLineColor(2);
	   pp[i+8]->SetStats(kFALSE);
	   pp[i+20]->SetLineColor(2);
	   pp[i+20]->SetStats(kFALSE);
	   pp[i+32]->SetLineColor(2);
	   pp[i+32]->SetStats(kFALSE);
	}
	for(i=6; i<12; ++i) {
	   pp[i+8]->SetLineColor(4);
	   pp[i+8]->SetStats(kFALSE);
	   pp[i+20]->SetLineColor(4);
	   pp[i+20]->SetStats(kFALSE);
	   pp[i+32]->SetLineColor(4);
	   pp[i+32]->SetStats(kFALSE);
	}
	   
	fscanf(ifp,"%f  %f  %f", &ysize, &xsize, &thick);
	zcen = thick/2.;
	printf("xsize/ysize/thick = %f/%f/%f \n", xsize, ysize, thick);
    fpix = false;
	if(thick > 286.) {fpix = true;}	   

	
    // Decide if this file corresponds to a CMSSW run 
    
	// Ask for external input into which pixels to join 
	
	printf("enter template ID \n");
	scanf("%d", &ID);
	printf("ID = %d \n", ID);
	
	
	// Initialize template store 
	
   std::vector< SiPixelTemplateStore > thePixelTemp_;
   SiPixelTemplate templ(thePixelTemp_);
	
	// Initialize template store, Pixelav 100V/300V simulation, +20C as thePixelTemp[6] 
   std::string templates_dir = "templates_dir/";
	
   templ.pushfile(ID,thePixelTemp_,templates_dir);
   templ.interpolate(ID, 0.f, 0.f, -1.f);
   qscale = templ.qscale();	
	
	// Initialize GenError store 
	
   std::vector< SiPixelGenErrorStore > thePixelGenErr_;
   SiPixelGenError gtempl(thePixelGenErr_);
	
   gtempl.pushfile(ID,thePixelGenErr_,templates_dir);	
   
	
	// Ask for speed info
	
	printf("enter algorithm speed (-2->5)\n");
	scanf("%d", &speed);
	
	printf("speed = %d\n", speed);
	
	// Ask for minimum probability
	
	printf("enter minimum probability\n");
	scanf("%lf", &probmin);
	
	printf("minimum probability = %lf\n", probmin);
		
		
	nevent=0;
	nbad = 0;
	tote = 0.;
	bade = 0.;
	locBz = -1.;
	
	static vector<float> sx(5,0.), sx2(5,0.), scx(5,0.), sy(5,0.), sy2(5,0.), scy(5,0.); 
	
	static vector<float> sxp(5,0.), sxp2(5,0.), scxp(5,0.), syp(5,0.), syp2(5,0.), scyp(5,0.); 
	
	static vector<float> sxc(5,0.), scxc(5,0.), sxc2(5,0.), syc(5,0.), scyc(5,0.), syc2(5,0.), nt(12,0.), nb(12,0.); 
	
	std::vector<std::pair<int, int> > pixlst;
	   
	
// Ask for external input into which pixels to join 

//	printf("enter the first y- and x-columns to double (-1 -1 = none)\n");
//	scanf("%d %d", &iyd, &ixd);
	iyd = -1; ixd = -1;
	
	printf("y/x double = %d/%d\n", iyd,ixd);
	
	int totale=0;  int goodx=0; int goody=0;
    
    ndcol = TYSIZE/2 + 1;
	std::vector<int> ndhit(ndcol, 0);
	
// Define tree quantities	
	
	int run_ = nfile;
	ULong64_t event_;
	int module_ = 7;
	int ladder_ = 6;
	int layer_ = 1;
	int isflipped_ = 1;
	float pt_;
	float p_;
	float eta_;
	float phi_;
	double chi2_;
	double ndof_;
	const int maxpix = 100;
	
	Pixinfo pixinfo_;
	
	Hit simhit_, trackhit_;
	
	Clust clust_;
	
	Rechit rechit_;
	
	TTRechit ttrechit_;
	
	//  Open a file to write the TTree structure
	
//  Create an input filename for this run 
	
	sprintf(ntupfile,"make_ntuple_output/pixelav_events_n%5.5d.root",nfile);
	TFile *f = new TFile(ntupfile,"RECREATE");
	f->cd();
    TTree *tree = new TTree("SiPixelLorentzAngleTree_", "pixelav simulated clusters");   
	tree->Branch("run", &run_,"run/I");	
	tree->Branch("event", &event_,"event/l");
	tree->Branch("module", &module_,"module/I");
	tree->Branch("ladder", &ladder_,"ladder/I");
	tree->Branch("layer", &layer_,"layer/I");
	tree->Branch("isflipped", &isflipped_,"isflipped/I");
	tree->Branch("pt", &pt_,"pt/F");
	tree->Branch("p", &p_,"p/F");
	tree->Branch("eta", &eta_,"eta/F");
	tree->Branch("phi", &phi_,"rphi/F");
	tree->Branch("chi2", &chi2_,"chi2/D");
	tree->Branch("ndof", &ndof_,"ndof/D");
	tree->Branch("trackhit", &trackhit_.x,"x/F:y/F:alpha/D:beta/D:gamma/D");
	tree->Branch("simhit", &simhit_.x,"x/F:y/F:alpha/D:beta/D:gamma/D");
	tree->Branch("npix", &pixinfo_.npix,"npix/I");
	tree->Branch("rowpix", pixinfo_.row,"rowpix[npix]/F");
	tree->Branch("colpix", pixinfo_.col,"colpix[npix]/F");
	tree->Branch("adc", pixinfo_.adc,"adc[npix]/F");
	tree->Branch("xpix", pixinfo_.x,"xpix[npix]/F");
	tree->Branch("ypix", pixinfo_.y,"ypix[npix]/F");
	tree->Branch("clust", &clust_.x,"x/F:y/F:charge/F:size_x/I:size_y/I:maxPixelCol/I:maxPixelRow/I:minPixelCol/I:minPixelRow/I");
	tree->Branch("rechit", &rechit_.x,"x/F:y/F");
	tree->Branch("ttrechit", &ttrechit_.x,"x/F:y/F:sigmax/F:sigmay/F:probx/F:proby/F:qbin/I");
	
/*  Determine current time */

	gettimeofday(&now0, &timz);
	   
// Loop until end of input file 

    while(fscanf(ifp,"%f %f %f %f %f %f %d", &pvec[0], &pvec[1], &pvec[2], &pvec[3], &pvec[4], &pvec[5], &neh) != EOF) {

// read the input cluster 

	   for (k=0; k < TXSIZE; ++k) {
		  fscanf(ifp,
			"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", 
			&pixin[k][0],&pixin[k][1],&pixin[k][2],&pixin[k][3],&pixin[k][4],&pixin[k][5],&pixin[k][6],&pixin[k][7],&pixin[k][8],&pixin[k][9],
			&pixin[k][10],&pixin[k][11],&pixin[k][12],&pixin[k][13],&pixin[k][14],&pixin[k][15],&pixin[k][16],&pixin[k][17],&pixin[k][18], 
			&pixin[k][19],&pixin[k][20]);
	   }
	   ++nevent;
       event_ = nevent;
       isflipped_ = 1;
       if(pvec[5] < 0.) {isflipped_ = 0;}
	   
// Add noise and analog response to cluster, reformat for flipped barrel coordinate system 
		 
       qmax = 0.;
       triplg(vgauss);
       pixlst.resize(0);
       for(i=0; i<ndcol; ++i) {ndhit[i] = 0;}
       icol = 0;
       if(vgauss[1] < 0.) {icol = 1;}
       for(j=0; j<TXSIZE; ++j) {
            triplg(wgauss);
            triplg(xgauss);
            triplg(ygauss);
            triplg(zgauss);
            for(i=0; i<TYSIZE; ++i) {
                bclust[j][i] = false;
                qin = (10.*pixin[j][i] + xgauss[i]*noise);
                rclust[TXSIZE-1-j][TYSIZE-1-i] = qin;
                if(qin < q100*(1.+wgauss[i]*q100_frac)) {
                    clust[TXSIZE-1-j][TYSIZE-1-i] = 0.;
                } else {
                    idcol = (TYSIZE-1-i+icol)/2;
                    ++ndhit[idcol];
                    if(non_linear == 0) {
                        qin *= (1.+gain_frac*ygauss[i]);
                        signal = (qin + zgauss[i]*readout_noise)/qscale;
                    } else {
						 adc = (double)((int)(p3+p2*tanh(p0*(qin + vcaloffst)/(7.0*vcal) - p1)));
						 signal = ((float)((1.+gain_frac*ygauss[i])*(vcal*gain*(adc-ped))) - vcaloffst + zgauss[i]*readout_noise)/qscale;

                    }	 
                    clust[TXSIZE-1-j][TYSIZE-1-i] = (1.+vgauss[0]*common_frac)*signal;
                    if(clust[TXSIZE-1-j][TYSIZE-1-i] > qmax) {
                        qmax = clust[TXSIZE-1-j][TYSIZE-1-i];
                        max.first = TXSIZE-1-j; max.second = TYSIZE-1-i;
                    }
                }
            }
        }
		 if(qmax < qseed) continue;
        
// Simulate the second, higher threshold in single dcol hits
        
        for(j=0; j<TXSIZE; ++j) {
            for(i=0; i<TYSIZE; ++i) {
                if(clust[j][i] > 0.) {
                    idcol = (i+icol)/2;
                    if(ndhit[idcol] == 1) {
                        // Apply higher threshold on single hits in double columns
                        if(rclust[j][i] < q101*(1.+wgauss[i]*q100_frac)) {
                            clust[j][i] = 0.;
                        }
                    }
                }
            }
        }
        
// Simulate clustering around maximum signal (seed)
		 
		 pixlst.push_back(max);
		 j=max.first; i=max.second;
		 bclust[j][i] = true;
		 
		 std::vector<std::pair<int, int> >::const_iterator pixIter, pixEnd;
		 
	 rescn: pixIter = pixlst.begin();
		 pixEnd = pixlst.end();
		 numadd = 0;
		 for ( ; pixIter != pixEnd; ++pixIter ) {
			 jmin = pixIter->first-1; 
			 jmax = jmin+3;
			 if(jmin < 0) {jmin = 0;}
			 if(jmax > TXSIZE) {jmax = TXSIZE;}
			 imin = pixIter->second-1;
			 imax = imin+3;
			 if(imin < 0) {imin = 0;}
			 if(imax > TYSIZE) {imax = TYSIZE;}
			 for(j=jmin; j<jmax; ++j) {
				 for(i=imin; i<imax; ++i) {
					 if(clust[j][i] > 0.) {
						 if(!bclust[j][i]) {
							 bclust[j][i] = true;
							 pixel.first = j; pixel.second = i;
							 pixlst.push_back(pixel);
							 ++numadd;
						 }
					 }
				 }
			 }
		 }
		 if(numadd > 0) goto rescn;
		
		for(j=0; j<TXSIZE; ++j) {
			for(i=0; i<TYSIZE; ++i) {
				cluster[j][i] = 0.;
			}
		}
		pixIter = pixlst.begin();
		pixEnd = pixlst.end();
		qclust=0.;
		npix = 0;
		colmin = 416; colmax = 0; imin = 0; imax = 0;
		rowmin = 160; rowmax = 0; jmin = 0; jmax = 0;
		for ( ; pixIter != pixEnd; ++pixIter ) {
			j = pixIter->first; 
			i = pixIter->second;
			cluster[j][i] = clust[j][i];
			qclust += clust[j][i];
// Calculate row and column numbers (remember that we will merge the first two rows/columns);
			row = j+81;
			col = i+209;
			pixinfo_.row[npix] = row;
			pixinfo_.col[npix] = col;
			if(row < rowmin) {rowmin = row; jmin=j;}
			if(row > rowmax) {rowmax = row; jmax = j;}
			if(col < colmin) {colmin = col; imin = i;}
			if(col > colmax) {colmax = col; imax = i;}			
			pixinfo_.adc[npix] = clust[j][i];
			pixinfo_.x[npix] = 0.0250+0.0100*j;
			pixinfo_.y[npix] = 0.0375+0.0150*i;
			++npix;
			if(npix == 100) break;
		}
		pixinfo_.npix = npix;
		clust_.x = (pixinfo_.x[jmin]+pixinfo_.x[jmax])/2.;
		clust_.y = (pixinfo_.y[imin]+pixinfo_.y[imax])/2.;
		clust_.charge = qclust/1000.;
		clust_.size_x = jmax - jmin + 1;
		clust_.size_y = imax - imin + 1;
		clust_.minPixelRow = rowmin;
		clust_.maxPixelRow = rowmax;
		clust_.minPixelCol = colmin;
		clust_.maxPixelCol = colmax;
		
		// Calculate the hit coordinates in the flipped coordinate system 
		
		yhit = -(pvec[0] + (zcen-pvec[2])*pvec[3]/pvec[5]);
		xhit = -(pvec[1] + (zcen-pvec[2])*pvec[4]/pvec[5]);
		
		// Do the template analysis on the cluster 
		
		cotalpha = pvec[4]/pvec[5];
		alpha = atan2((double)pvec[5], (double)pvec[4]);
		cotbeta = pvec[3]/pvec[5];
		beta = atan2((double)pvec[5], (double)pvec[3]);
// Change this to agree with Mirena's tree-builder code
		gamma = atan2((double)pvec[4], (double)pvec[3]);
		etas = log((double)(-cotbeta+sqrt((double)(1.+cotbeta*cotbeta))));
		eta_ = etas;
		p_ = sqrt((double)(pvec[3]*pvec[3]+pvec[4]*pvec[4]+pvec[5]*pvec[5]));
		pt_ = p_ * fabsf(pvec[5])/sqrt(pvec[3]*pvec[3]+pvec[5]*pvec[5]);
		eta = fabs(etas);
		weight = 1./cosh((double)eta);
		tote += weight;
		etabin = (int)(eta/0.25);
		if(etabin > 11) {etabin = 11;}
		++nt[etabin];
		
// Work only on a single pixel part of the sensor (row/column 81/209 and beyond)
		
		for(i=0; i<TYSIZE; ++i) {
			ydouble[i] = false;
		}
		
		for(j=0; j<TXSIZE; ++j) {
			xdouble[j] = false;
		}
	   
// Do the template analysis on the cluster 

// Do the template analysis on the cluster 
       SiPixelTemplateReco::ClusMatrix clusterPayload{&cluster[0][0], xdouble, ydouble, mrow,mcol};
       SiPixelTemplateReco::ClusMatrix clusterPayloadC = clusterPayload;
	   locBx = 1.;
       if(cotbeta < 0.) locBx = -1.;
       locBz = locBx;
       if(cotalpha < 0.) locBz = -locBx;
       ierr = PixelTempReco1D(ID, cotalpha, cotbeta, locBz, locBx, clusterPayload, templ, yrec, sigmay, proby, xrec, sigmax, probx, qbin, speed, probQ);		 
	   if(ierr != 0) {
	      ++nbad; ++nb[etabin]; bade +=weight;
	      printf("reconstruction failed with error %d \n", ierr);
		} else {
	   
// Check resolution and weights 
           if(qbin > 3) {
			   qbin = 4;
		      //printf(" qbin = %d \n", qbin);
              ++nbin[qbin];
		   }
           ++nbin[qbin];
		   yhitcm = 0.0001*(yhit+((float)(TYSIZE)/2.+2.0)*ysize);
		   dy = yrec - (TYSIZE/2)*ysize - yhit;
		   yreccm = 0.0001*(yrec+2.5*ysize);
		   ++totale;
		   if(probx > 1.e-2) {++goodx;}else{hp[20]->Fill((double)cotbeta);hp[21]->Fill(alpha);}
		   if(proby > 1.e-2) {++goody;}
		   
		   log10proby = log10((double)proby); log10probx = log10((double)probx);
		   sy[qbin] += dy; sy2[qbin] += dy*dy; scy[qbin] += dy*dy/(sigmay*sigmay);
		   if(proby > 1.e-2) {syp[qbin] += dy; syp2[qbin] += dy*dy; scyp[qbin] += dy*dy/(sigmay*sigmay);}
		   hp[0]->Fill(dy);
		   hp[1+qbin]->Fill(dy);
		   hp[5]->Fill(log10proby);
		   hp[6+qbin]->Fill(log10proby);
	       if(qbin == 1) {pp[0]->Fill(eta, dy);}
		   if(qbin > 1) {pp[2]->Fill(eta, dy);}
		   proba = probx*proby;
		   probXY = proba*(1.-log(proba));
		   if(probXY < probmin) continue;			
		   xhitcm = 0.0001*(xhit+((float)(TXSIZE)/2.+2.)*xsize);
		   dx = xrec - (TXSIZE/2)*xsize - xhit;
		   xreccm = 0.0001*(xrec+2.5*xsize);
			
// Fill ntuple quantities
			simhit_.x = xhitcm;
			simhit_.y = yhitcm;
			simhit_.alpha = alpha;
			simhit_.beta = beta;
         simhit_.gamma = gamma;
         trackhit_.x = xreccm;
         trackhit_.y = yreccm;
         trackhit_.alpha = alpha;
         trackhit_.beta = beta;
         trackhit_.gamma = gamma;
			ttrechit_.x = xreccm;
			ttrechit_.y = yreccm;
			ttrechit_.sigmax = 0.0001*sigmax;
			ttrechit_.sigmay = 0.0001*sigmay;
			ttrechit_.probx = probx;
			ttrechit_.proby = proby;
			ttrechit_.qbin = qbin;
			phi_ = 0.8976 + xhitcm/11.;
			
			sx[qbin] += dx; sx2[qbin] += dx*dx; scx[qbin] += dx*dx/(sigmax*sigmax);
		   if(probx > 1.e-2) {sxp[qbin] += dx; sxp2[qbin] += dx*dx; scxp[qbin] += dx*dx/(sigmax*sigmax);}
		   hp[10]->Fill(dx);
		   hp[11+qbin]->Fill(dx);
		   hp[15]->Fill(log10probx);
		   hp[16+qbin]->Fill(log10probx);
		   if(qbin == 1) {pp[4]->Fill(eta, dx);}
		   if(qbin > 1) {pp[6]->Fill(eta, dx);}
 		   if(qbin == 1) {pp[44]->Fill(alpha, dx);}
		   if(qbin > 1) {pp[46]->Fill(alpha, dx);}
	        locBx = 1.;
            if(cotbeta < 0.) locBx = -1.;
            locBz = locBx;
            if(cotalpha < 0.) locBz = -locBx;            
            ierr = PixelGeneric2D(ID, cotalpha, cotbeta, locBz, locBx, clusterPayloadC, gtempl, ycmssw, sycmssw, xcmssw, sxcmssw, nypix, nxpix, yfrac, xfrac);
		   dyc = ycmssw - (TYSIZE/2)*ysize - yhit;
		   ycmsswcm = 0.0001*(ycmssw+2.5*ysize);
		   syc[qbin] += dyc; syc2[qbin] += dyc*dyc; scyc[qbin] += dyc*dyc/(sycmssw*sycmssw);
			dxc = xcmssw - (TXSIZE/2)*xsize - xhit;
			xcmsswcm = 0.0001*(xcmssw+2.5*xsize);
			rechit_.x = xcmsswcm;
			rechit_.y = ycmsswcm;
			chi2_ = (double)(dxc*dxc/(sxcmssw*sxcmssw) + dyc*dyc/(sycmssw*sycmssw));
			ndof_ = 2.;
			tree->Fill();
		   sxc[qbin] += dxc; sxc2[qbin] += dxc*dxc; scxc[qbin] += dxc*dxc/(sxcmssw*sxcmssw);
		   hp[22]->Fill(dyc);
		   hp[23]->Fill(dxc);
		   if(qbin == 1) {pp[1]->Fill(eta, dyc);}
		   if(qbin > 1) {pp[3]->Fill(eta, dyc);}
		   if(qbin == 1) {pp[5]->Fill(eta, dxc);}
		   if(qbin > 1) {pp[7]->Fill(eta, dxc);}
		   if(qbin == 1) {pp[45]->Fill(alpha, dxc);}
		   if(qbin > 1) {pp[47]->Fill(alpha, dxc);}
		   if(nypix == 2 && yfrac >= 0.) {pp[8]->Fill(yfrac, dyc);}
		   if(nypix == 3 && yfrac >= 0.) {pp[9]->Fill(yfrac, dyc);}
		   if(nypix == 4 && yfrac >= 0.) {pp[10]->Fill(yfrac, dyc);}
		   if(nxpix == 2 && xfrac >= 0.) {pp[11]->Fill(xfrac, dxc);}
		   if(nxpix == 3 && xfrac >= 0.) {pp[12]->Fill(xfrac, dxc);}
		   if(nxpix == 4 && xfrac >= 0.) {pp[13]->Fill(xfrac, dxc);}
		   if(nypix == 2 && yfrac >= 0.) {pp[14]->Fill(yfrac, dy);}
		   if(nypix == 3 && yfrac >= 0.) {pp[15]->Fill(yfrac, dy);}
		   if(nypix == 4 && yfrac >= 0.) {pp[16]->Fill(yfrac, dy);}
		   if(nxpix == 2 && xfrac >= 0.) {pp[17]->Fill(xfrac, dx);}
		   if(nxpix == 3 && xfrac >= 0.) {pp[18]->Fill(xfrac, dx);}
		   if(nxpix == 4 && xfrac >= 0.) {pp[19]->Fill(xfrac, dx);}
		   if(qbin==1) {
		      if(nypix == 2 && yfrac >= 0.) {pp[20]->Fill(yfrac, dyc);}
		      if(nypix == 3 && yfrac >= 0.) {pp[21]->Fill(yfrac, dyc);}
		      if(nypix == 4 && yfrac >= 0.) {pp[22]->Fill(yfrac, dyc);}
		      if(nxpix == 2 && xfrac >= 0.) {pp[23]->Fill(xfrac, dxc);}
		      if(nxpix == 3 && xfrac >= 0.) {pp[24]->Fill(xfrac, dxc);}
		      if(nxpix == 4 && xfrac >= 0.) {pp[25]->Fill(xfrac, dxc);}
		      if(nypix == 2 && yfrac >= 0.) {pp[26]->Fill(yfrac, dy);}
		      if(nypix == 3 && yfrac >= 0.) {pp[27]->Fill(yfrac, dy);}
		      if(nypix == 4 && yfrac >= 0.) {pp[28]->Fill(yfrac, dy);}
		      if(nxpix == 2 && xfrac >= 0.) {pp[29]->Fill(xfrac, dx);}
		      if(nxpix == 3 && xfrac >= 0.) {pp[30]->Fill(xfrac, dx);}
		      if(nxpix == 4 && xfrac >= 0.) {pp[31]->Fill(xfrac, dx);}
           }		   
		   if(qbin > 1) {
		      if(nypix == 2 && yfrac >= 0.) {pp[32]->Fill(yfrac, dyc);}
		      if(nypix == 3 && yfrac >= 0.) {pp[33]->Fill(yfrac, dyc);}
		      if(nypix == 4 && yfrac >= 0.) {pp[34]->Fill(yfrac, dyc);}
		      if(nxpix == 2 && xfrac >= 0.) {pp[35]->Fill(xfrac, dxc);}
		      if(nxpix == 3 && xfrac >= 0.) {pp[36]->Fill(xfrac, dxc);}
		      if(nxpix == 4 && xfrac >= 0.) {pp[37]->Fill(xfrac, dxc);}
		      if(nypix == 2 && yfrac >= 0.) {pp[38]->Fill(yfrac, dy);}
		      if(nypix == 3 && yfrac >= 0.) {pp[39]->Fill(yfrac, dy);}
		      if(nypix == 4 && yfrac >= 0.) {pp[40]->Fill(yfrac, dy);}
		      if(nxpix == 2 && xfrac >= 0.) {pp[41]->Fill(xfrac, dx);}
		      if(nxpix == 3 && xfrac >= 0.) {pp[42]->Fill(xfrac, dx);}
		      if(nxpix == 4 && xfrac >= 0.) {pp[43]->Fill(xfrac, dx);}
           }		   
	    }
		
   }
	
	tree->Write();
	tree->Show(10);
	tree->Show(50);
	f->Close();
   
/*  Determine current time */

   gettimeofday(&now1, &timz);
   deltas = now1.tv_sec - now0.tv_sec;
   deltaus = now1.tv_usec - now0.tv_usec;
   deltat = ((double)deltaus)/1000000.;
   deltat += (double)deltas;
   printf("ellapsed time = %f seconds \n", deltat);
   
   printf(" total events = %d, probx > 10^{-3} = %d, proby > 10^{-3} = %d \n", totale, goodx, goody);
   printf(" low q failures = %d, malformed clusters = %d \n", nbin[4], nbad);
	   
   for(j=0; j<5; ++j) {
      sy[j] /= (float)nbin[j]; sy2[j] /= (float)nbin[j]; scy[j] /= (float)nbin[j];
      sy2[j] = sqrt((double)(sy2[j] - sy[j]*sy[j]));
      printf(" avg y residual[%1d] = %f +- %f, avg y chi^2 = %f \n", j, sy[j], sy2[j], scy[j]);       
      sx[j] /= (float)nbin[j]; sx2[j] /= (float)nbin[j]; scx[j] /= (float)nbin[j];
      sx2[j] = sqrt((double)(sx2[j] - sx[j]*sx[j]));
      printf(" avg x residual[%1d] = %f +- %f, avg x chi^2 = %f \n", j, sx[j], sx2[j], scx[j]); 
  }
   printf(" After 10^{-2} probability cuts: \n");
	   
   for(j=0; j<5; ++j) {
      syp[j] /= (float)nbin[j]; syp2[j] /= (float)nbin[j]; scyp[j] /= (float)nbin[j];
      syp2[j] = sqrt((double)(syp2[j] - syp[j]*syp[j]));
      printf(" avg y residual[%1d] = %f +- %f, avg y chi^2 = %f \n", j, syp[j], syp2[j], scyp[j]);       
      sxp[j] /= (float)nbin[j]; sxp2[j] /= (float)nbin[j]; scxp[j] /= (float)nbin[j];
      sxp2[j] = sqrt((double)(sxp2[j] - sxp[j]*sxp[j]));
      printf(" avg x residual[%1d] = %f +- %f, avg x chi^2 = %f \n", j, sxp[j], sxp2[j], scxp[j]); 
  }
  
   printf(" CMSSW algorithm \n");
	   
   for(j=0; j<5; ++j) {
      syc[j] /= (float)nbin[j]; syc2[j] /= (float)nbin[j]; scyc[j] /= (float)nbin[j];
      syc2[j] = sqrt((double)(syc2[j] - syc[j]*syc[j]));
      printf(" avg y residual[%1d] = %f +- %f, avg y chi^2 = %f \n", j, syc[j], syc2[j], scyc[j]);       
      sxc[j] /= (float)nbin[j]; sxc2[j] /= (float)nbin[j]; scxc[j] /= (float)nbin[j];
      sxc2[j] = sqrt((double)(sxc2[j] - sxc[j]*sxc[j]));
      printf(" avg x residual[%1d] = %f +- %f, avg x chi^2 = %f \n", j, sxc[j], sxc2[j], scxc[j]); 
  }
  
// Make resolution graphs from profile errors

  float binx, biny, bino;
  printf("\n Fraction malformed clusters = %lf\n", bade/tote);
  for(j=0; j<11; ++j) {
	binx=0.125+j*0.25; biny=nb[j]/nt[j];
	printf("%f  %f \n",binx, biny);
  }
  for(i=0; i<8; ++i) {
  printf("\n Profile errors %d, %s \n",i, pp[i]->GetTitle());
     for(j=0; j<11; ++j) {
	    binx=pp[i]->GetBinCenter(j+1); biny=pp[i]->GetBinError(j+1);
        printf("%f  %f \n",binx, biny);
	 }
  }
  
   for(i=44; i<48; ++i) {
  printf("\n Profile errors %d, %s \n",i, pp[i]->GetTitle());
     for(j=0; j<25; ++j) {
	    binx=pp[i]->GetBinCenter(j+1); biny=pp[i]->GetBinError(j+1), bino=pp[i]->GetBinContent(j+1);
        printf("%f  %f %f \n",binx, biny, bino);
	 }
  }
 
/*
 * Histograms plotting
 */
	TCanvas c2("c2", header);
	c2.SetFillStyle(4000);
   for(i=0; i<5; ++i) {hp[i]->Fit("gaus"); hp[i+10]->Fit("gaus");}
	for(i=22; i<24; ++i) {hp[i]->Fit("gaus");}
	
//  Create an output filename for this run 

   sprintf(outfile0,"make_ntuple_output/pixel_histos%5.5d.pdf[",nfile);
   sprintf(outfile1,"make_ntuple_output/pixel_histos%5.5d.pdf",nfile);
   sprintf(outfile2,"make_ntuple_output/pixel_histos%5.5d.pdf]",nfile);
   c2.Print(outfile0);
   for(i=0; i<24; ++i) {
     hp[i]->Draw();
     c2.Print(outfile1);
   }
   for(i=0; i<8; ++i) {
     pp[i]->Draw();
     c2.Print(outfile1);
   }
   for(i=8; i<14; ++i) {
     pp[i]->Draw();
	 pp[i+6]->Draw("same");
     c2.Print(outfile1);
   }
   for(i=20; i<26; ++i) {
     pp[i]->Draw();
	 pp[i+6]->Draw("same");
     c2.Print(outfile1);
   }
   for(i=32; i<38; ++i) {
     pp[i]->Draw();
	 pp[i+6]->Draw("same");
     c2.Print(outfile1);
   }
   c2.Print(outfile2);
   
return 0;
} // MAIN__ 




// ***************************************************************** 
//! Calculate 21 gaussianly-distributed random numbers.
//! \param x - a vector holding 21 random numbers generated with rms = 1.            
// ***************************************************************** 
int triplg(std::vector<float>& x)
{
    // Initialized data 

    static int fcall = -1;

    // Local variables 
    static float r1, r2;
    static int i__;
    static float r__;
    static int ibase;
    static std::vector<float> rbuff(TYTEN);
    static float twopi;
    double arg, phi;



    // Function Body 

//  Initalize the parameters 

    if (fcall) {
	twopi = 2.*acos((double)-1.);
	ibase = TYTEN;
	fcall = 0;
    }

//  If all random numbers used up, generate 210 more 

    if (ibase == TYTEN) {
	   for (i__ = 0; i__ < TYTEN-1; i__ += 2) {
	      r1 = ((float)random())/((float)RAND_MAX);
	      r2 = ((float)random())/((float)RAND_MAX);
	      arg = (double)(1. - r1);
	      if (arg < 1.e-30) {arg = 1.e-30;}
	      r__ = sqrt(log(arg) * (-2.));
	      phi = twopi * r2;
	      rbuff[i__] = r__ * cos(phi);
	      rbuff[i__+1] = r__ * sin(phi);
	   }
	   ibase = 0;
    }
    for (i__ = 0; i__ < TYSIZE; ++i__) {
	   x[i__] = rbuff[ibase + i__];
    }
    ibase += TYSIZE;
    return 0;
} // triplg 

