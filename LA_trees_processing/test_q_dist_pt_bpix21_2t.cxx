//! \file template_code6d.cpp
//!
//! Apply standard CMSSW pixel reconstruction to pixelav hits 
//! Compare with Lorentz data clusters
//! Try applying probability cut first
//! Parallel stream to try larger readout threshold
//! Use consistent definition of scale factor [q_correct = q_cmssw*qscale]
//! Change smearing to pixel-by-pixel and common-mode
//! Read new pt summary format
//! Add second threshold to simulation
//! Rearrange the histograms
//! Add layer info 
//! Add Plus/Minus end of layers
//! Add flipped/unflipped modules


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
static int theVerboseLevel = {2};
#include "SiPixelTemplateReco.cc"
#include "VVIObj.cc"
//#include "PixelGeneric2D.cc"

using namespace std;

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"

#define TYSIZE 21
#define TYTEN 210 // = 10*TYSIZE
#define BYSIZE TYSIZE+4
#define BHY 12 // = BYSIZE/2
#define BYM1 TYSIZE+3
#define BYM2 TYSIZE+2
#define BYM3 TYSIZE+1
#define TXSIZE 13
#define BXSIZE TXSIZE+4
#define BHX 8 // = BXSIZE/2
#define BXM1 TXSIZE+3
#define BXM2 TXSIZE+2
#define BXM3 TXSIZE+1


// Global definitions 


// Main program  

int main(int argc, char *argv[])
{
   
    std::vector<float> pvec(6), wgauss(TYSIZE), vgauss(TYSIZE), xgauss(TYSIZE), ygauss(TYSIZE), zgauss(TYSIZE);
    static bool fpix;
    float pixin[TXSIZE][TYSIZE];
    bool ydouble[TYSIZE], xdouble[TXSIZE];
        std::vector<int> nd1(4,0), nd2(4,0), nd3(4,0), nd4(4,0);
    static float thick, xsize, ysize, noise, zcen, gain_frac, readout_noise, q100_frac, common_frac;
	static float xhit, yhit, xrec, yrec, sigmax, sigmay, probx, proby, signal, cotalpha, cotbeta, probQH, probQ, tkp, tkpmin, qscale, locBx, locBz;  
    static int sizex, sizey, ntrack, ntrack2, ngood, ngood2, nbad, icol, ndcol;
    static int ndata, nfile, neh, nevent, ID, non_linear, layer, ID0, end, ltpt; 
    int i, j, k, ierr, qbin, numadd, flipped, module, idcol, lay, flp;
    int mrow = TXSIZE, mcol = TYSIZE;
	double dx, dy, eta, proba, probb, qtotal, pcut, log10px, log10py, log10pQ, qnorm, adc;
	static int iyd, ixd, speed, qbcut, imax, imin, jmax, jmin;	
	static float q100, q101, q50, q10, qmax; 
	const double gain = 3.19;
	const double ped = 16.46;
	const double p0 = 0.01218;
	const double p1 = 0.711;
	const double p2 = 203.;
	const double p3 = 148.;	
	static double vcal = 47.;	
	static double vcaloffst = 60.;
	float sigtmp, qin;
    static char infile[80], header[80], title[200], c, outfile0[80], outfile1[80], outfile2[80];
	int triplg(std::vector<float>&);
	int triplu(std::vector<float>&);
//	int random(void);
	float cluster[TXSIZE][TYSIZE], clust[TXSIZE][TYSIZE], rclust[TXSIZE][TYSIZE];
	bool bclust[TXSIZE][TYSIZE];
    std::pair<int, int> pixel, max;

    FILE *ifp;
	
   struct timeval now0, now1;
   struct timezone timz;
   long deltas, deltaus;
   double deltat;
        			     
	
	// Ask for external input into which pixels to join 
	
	printf("enter probability cut, qbin cut (4/5), momentum cut (GeV), Template ID, layer (1-6,0=all,7=all new), Min/Pls/Both (0/1/2), unf/flp/both (0/1/2), lt/pt (0/1) \n");
	scanf("%lf %d %f %d %d %d %d %d", &pcut, &qbcut, &tkpmin, &ID0, &layer, &end, &flp, &ltpt);
//	printf("probability cut = %lf, qbin cut = %d, momentum cut = %f, ID0 = %d, layer = %d, end = %d, flipped = %d, lt/pt = %d \n", pcut, qbcut, tkpmin, ID0, layer, end, flp, ltpt);
	
    //	double  halfxs=300.;
	int nx=120;	
	gStyle->SetOptFit(101);
	gStyle->SetHistLineWidth(2);
	static vector<TH1F*> hp(49);
    sprintf(title,"Detector Cluster Charge (all clust, p*100>%6.4lf ); Cluster Charge (e)", pcut*100.);
    hp[0] = new TH1F("h701",title,nx,0.,300000.);
    sprintf(title,"Detector Cluster Charge (size>1, 0<qbin<%d, p*100>%6.4lf ); Cluster Charge (e)", qbcut, pcut*100.);
    hp[2] = new TH1F("h703",title,nx,0.,300000.);
    sprintf(title,"Detector qbin (all clust, p*100>%6.4lf ); Qbin", pcut*100.);
    hp[4] = new TH1F("h101",title,6,-0.5,5.5);
    sprintf(title,"Detector qbin (size>1, p*100>%6.4lf ); Qbin", pcut*100.);
    hp[6] = new TH1F("h103",title,6,-0.5,5.5);
    sprintf(title,"Detector Px (All Q, p*100>%6.4lf ); Prob_x", pcut*100.);
    hp[8] = new TH1F("h105",title,100,0.0,1.0);
    sprintf(title,"Detector Px (Qbin>0, p*100>%6.4lf ); Prob_x", pcut*100.);
    hp[10] = new TH1F("h107",title,100,0.0,1.0);
    sprintf(title,"Detector Py (All Q, p*100>%6.4lf ); Prob_y", pcut*100.);
    hp[12] = new TH1F("h109",title,100,0.0,1.0);
    sprintf(title,"Detector Py (Qbin>0, p*100>%6.4lf ); Prob_y", pcut*100.);
    hp[14] = new TH1F("h111",title,100,0.0,1.0);
    sprintf(title,"Detector log10(Px) (All Q, p*100>%6.4lf ); log10(P_x)", pcut*100.);
    hp[16] = new TH1F("h113",title,120,-6.0,0.0);
    sprintf(title,"Detector log10(Px) (Qbin>0, p*100>%6.4lf ); log10(P_x)", pcut*100.);
    hp[18] = new TH1F("h115",title,120,-6.0,0.0);
    sprintf(title,"Detector log10(Py) (All Q, p*100>%6.4lf ); log10(P_y)", pcut*100.);
    hp[20] = new TH1F("h117",title,120,-6.0,0.0);
    sprintf(title,"Detector log10(Py) (Qbin>0, p*100>%6.4lf ); log10(P_y)", pcut*100.);
    hp[22] = new TH1F("h119",title,120,-6.0,0.0);
    sprintf(title,"Detector Cluster Charge (size>1, qbin<%d, p*100>%6.4lf ); Cluster Charge (e)", qbcut, pcut*100.);
    hp[24] = new TH1F("h705",title,nx,0.,300000.);
    sprintf(title,"Detector Cluster Charge (size=1, qbin<%d, p*100>%6.4lf ); Cluster Charge (e)", qbcut, pcut*100.);
    hp[26] = new TH1F("h707",title,75,0.,75000.);
    sprintf(title,"Normalized Detector Cluster Charge (size>1, qbin<%d, p*100>%6.4lf ); Normalized Cluster Charge (e)", qbcut, pcut*100.);
    hp[28] = new TH1F("h709",title,nx,0.,80000.);
    sprintf(title,"Detector PQ (All Q, p*100>%6.4lf ); Prob_Q", pcut*100.);
    hp[30] = new TH1F("h121",title,100,0.0,1.0);
    sprintf(title,"Detector PQ (Qbin>0, p*100>%6.4lf ); Prob_Q", pcut*100.);
    hp[32] = new TH1F("h123",title,100,0.0,1.0);
    sprintf(title,"Detector log10(PQ) (All Q, p*100>%6.4lf ); log10(P_Q)", pcut*100.);
    hp[34] = new TH1F("h125",title,100,-5.0,0.0);
    sprintf(title,"Detector log10(PQ) (Qbin>0, p*100>%6.4lf ); log10(P_Q)", pcut*100.);
    hp[36] = new TH1F("h127",title,100,-5.0,0.0);
    sprintf(title,"Detector sizex (p*100>%6.4lf ); sizex (pix)", pcut*100.);
    hp[38] = new TH1F("h129",title,10, 0.5,10.5);
    sprintf(title,"Detector sizey (p*100>%6.4lf ); sizey (pix)", pcut*100.);
    hp[40] = new TH1F("h131",title,20,0.5,20.5);
    sprintf(title,"Normalized Detector Cluster Charge (all sizes, p*100>%6.4lf ); Normalized Cluster Charge (e)", pcut*100.);
    hp[42] = new TH1F("h711",title,nx,0.,80000.);
    sprintf(title,"Detector sizex (qbin>0, p*100>%6.4lf ); sizex (pix)", pcut*100.);
    hp[44] = new TH1F("h133",title,10, 0.5,10.5);
    sprintf(title,"Detector sizey (qbin>0, p*100>%6.4lf ); sizey (pix)", pcut*100.);
    hp[46] = new TH1F("h135",title,20,0.5,20.5);
	
	
	// Set style for the the histograms	
	
	hp[0]->SetLineColor(2);
	//	   hp[0]->SetFillColor(kBlue);
	hp[2]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[4]->SetLineColor(2);
	//	   hp[0]->SetFillColor(kBlue);
	hp[6]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[8]->SetLineColor(2);
	//	   hp[0]->SetFillColor(kBlue);
	hp[10]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[12]->SetLineColor(2);
	//	   hp[0]->SetFillColor(kBlue);
	hp[14]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[16]->SetLineColor(2);
	//	   hp[0]->SetFillColor(kBlue);
	hp[18]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[20]->SetLineColor(2);
	//	   hp[0]->SetFillColor(kBlue);
	hp[22]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[24]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[26]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[28]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[30]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[32]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[34]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[36]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[38]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[40]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[42]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[44]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	hp[46]->SetLineColor(2);
	//	   hp[2]->SetFillColor(kBlue);
	
	
	
	
//  Read which data and inputs to use (use c file i/o which is much faster than c++ i/o) 
	
	ifp = fopen("q_dist_2t.txt", "r");
	if (ifp==NULL) {
      printf("no pixel initialization file found/n");
      return 0;
	}
	
	fscanf(ifp,"%d %d %f %f %f %f %f %f %f %d", &ndata, &nfile, &noise, &q100, &q101, &q100_frac, &common_frac, &gain_frac, &readout_noise, &non_linear);
	fclose(ifp);
//	printf("data file %d, mc file %d noise = %f, threshold0 = %f, threshold1 = %f, rms threshold frac = %f, common_frac = %f, gain fraction = %f, readout noise = %f, nonlinear_resp = %d \n", ndata, nfile, noise, q100, q101, q100_frac, common_frac, gain_frac, readout_noise, non_linear);
	
//  Calculate 50% of threshold in q units and enc noise in adc units
	
	q50=0.5*q100;
	q10=0.2*q50;
	
//  Create an input data file name for this run 
	
	if(ltpt == 0) {sprintf(infile,"Efield_output/lt_bpix%6.6d.txt",ndata);} else {sprintf(infile,"pt_bpix%6.6d.txt",ndata);}
	
	//  Open input file and read header info 
	
	ifp = fopen(infile, "r");
	if (ifp==NULL) {
		printf("no pixel data file found/n");
		return 0;
	}
	
	ntrack = 0; ntrack2 = 0;
    while(fscanf(ifp,"%f %f %lf %d %d %d %f %f %lf %lf %f %f %d %d %d", &cotalpha, &cotbeta, &qtotal, &sizex, &sizey, &qbin, &probx, &proby, &dx, &dy, &probQH, &tkp, &flipped, &module, &lay) != EOF) {
//        printf("%f %f %lf %d %d %d %f %f %lf %lf \n", cotalpha, cotbeta, qtotal, sizex, sizey, qbin, probx, proby, dx, dy);
		 if(probQH > 0.5) {
			 probQ = 2.*(1. - probQH);
		 } else {
			 probQ = 2.*probQH;
		 }
       if(layer == 7 && lay > 3) lay = 7;
		 if(tkp < tkpmin) continue;
        if(layer > 0 && lay != layer) continue;
        if(end == 0 && module > 0) continue;
        if(end == 1 && module < 0) continue;
        if(flp < 2 && flp != flipped) continue;
		++nd1[0];
		if(sizex > 1 || sizey > 1) {++nd1[1];}
		if(probx < 1.e-6) probx = 1.e-6;
		if(proby < 1.e-6) proby = 1.e-6;
		proba = (double)(probx*proby);
		probb = proba*(1.-log(proba));
		log10px=log10((double)probx);
		log10py=log10((double)proby);
		log10pQ=log10((double)probQ);
		qnorm = qtotal/sqrt((double)(1.+cotbeta*cotbeta+cotalpha*cotalpha));			   
		hp[42]->Fill(qnorm);
		if(qbin < qbcut) {
		   ++nd2[0];
		   if(sizex > 1 || sizey > 1) {++nd2[1];}
		}
		if(probb > pcut) {
		   ++ntrack; 
//		   printf("qtotal = %lf \n", qtotal);
		   hp[0]->Fill(qtotal);
		   hp[4]->Fill((double)qbin);
			hp[8]->Fill((double)probx);
			hp[12]->Fill((double)proby);
			hp[30]->Fill((double)probQ);
			hp[16]->Fill(log10px);
			hp[20]->Fill(log10py);
			hp[34]->Fill(log10pQ);
		   if(sizex == 1 && sizey == 1) {
			   if(qbin<qbcut) {hp[26]->Fill(qtotal);}
		   }
		   ++nd3[0];
			hp[38]->Fill((double)sizex);
			hp[40]->Fill((double)sizey);
			if(qbin > 0) {
				hp[44]->Fill((double)sizex);
				hp[46]->Fill((double)sizey);
			}
		   if(qbin > 0) {
			   hp[10]->Fill((double)probx);
			   hp[14]->Fill((double)proby);
			   hp[32]->Fill((double)probQ);
			   hp[18]->Fill(log10px);
			   hp[22]->Fill(log10py);
			   hp[36]->Fill(log10pQ);
			}
		   if(sizex > 1 || sizey > 1) {
			  ++nd3[1];
			  hp[6]->Fill((double)qbin);			
			   qnorm = qtotal/sqrt((double)(1.+cotbeta*cotbeta+cotalpha*cotalpha));			   
			  if(qbin < qbcut) {hp[28]->Fill(qnorm);}
			  if(qbin < qbcut && qbin > 0) {hp[2]->Fill(qtotal); ++ntrack2;}
			  if(qbin<qbcut) {hp[24]->Fill(qtotal);}
		  }
		  if(qbin < qbcut) {
			  ++nd4[0];
			  if(sizex > 1 || sizey > 1) {++nd4[1];}
		  }
	   }
	}
    fclose(ifp);
    printf("mc file number = %d \n", nfile);
    
	printf("number of clusters = %d, number of n>1 clusters = %d \n", ntrack, ntrack2);
	
	
	//  Create an input filename for this run 
	
    if(nfile < 10000) {
        
        sprintf(infile,"simulated_clusters/template_events_d%4.4d.out",nfile);
        
    } else {
        
        sprintf(infile,"simulated_clusters/template_events_d%5.5d.out",nfile);
        
    }
	
	//  Open input file and read header info 
	
	ifp = fopen(infile, "r");
	if (ifp==NULL) {
      printf("no pixel simuated data file found/n");
      return 0;
	}
	
// Read-in a header string first and print it    
    
    for (i=0; (c=getc(ifp)) != '\n'; ++i) {
       if(i < 79) {header[i] = c;}
    }
	if(i > 78) {i=78;}
	header[i+1] ='\0';
    printf("%s\n", header);
	
	
// Make a 2D histogram
//	static vector<TH2F*> h2(1);
//	h2[0] = new TH2F("h2101","cot(#alpha); cot(#beta)",40,-4.0,4.0,40,-2.0,2.0);	
	   
	fscanf(ifp,"%f  %f  %f", &ysize, &xsize, &thick);
	zcen = thick/2.;
	printf("xsize/ysize/thick = %f/%f/%f \n", xsize, ysize, thick);
    fpix = false;
	if(thick > 286.) {fpix = true;}
	
	nevent=0;
		
	std::vector<std::pair<int, int> > pixlst;
	
	
	// Decide if this file corresponds to a CMSSW run 
	
	if(layer == 1 || layer == 5) {
	   vcal = 50.;	
	   vcaloffst = 670;
	} 
    
    ID = ID0;
	
	if(fpix) {++ID; printf("ID = %d, fpix \n", ID);} else {printf("ID = %d, barrel, vcal = %lf, offset = %lf \n", ID, vcal, vcaloffst);}
	
	// Initialize template store 
	
   std::vector< SiPixelTemplateStore > thePixelTemp_;
   SiPixelTemplate templ(thePixelTemp_);
	
	// Initialize template store, Pixelav 100V/300V simulation, +20C as thePixelTemp[6] 
   std::string templates_dir = "templates_dir/";
   templ.pushfile(ID,thePixelTemp_,templates_dir);	
   templ.interpolate(ID, 0.f, 0.f, -1.f);
   qscale = templ.qscale();

	
	sprintf(title,"Pixelav Cluster Charge (all clust, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Cluster Charge (e)", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[1] = new TH1F("h702",title,nx,0.,300000.);      
	sprintf(title,"Pixelav Cluster Charge (size>1, 0<qbin<%d, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Cluster Charge (e)", qbcut, pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[3] = new TH1F("h704",title,nx,0.,300000.);      
	sprintf(title,"Pixelav qbin (all clust, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Qbin", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[5] = new TH1F("h102",title,6,-0.5,5.5);
	sprintf(title,"Pixelav qbin (size>1, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Qbin", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[7] = new TH1F("h104",title,6,-0.5,5.5);
	sprintf(title,"Pixelav Px (All Q, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Prob_x", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[9] = new TH1F("h106",title,100,0.0,1.0);
	sprintf(title,"Pixelav Px (Qbin>0, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Prob_x", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[11] = new TH1F("h108",title,100,0.0,1.0);
	sprintf(title,"Pixelav Py (All Q, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Prob_y", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[13] = new TH1F("h110",title,100,0.0,1.0);
	sprintf(title,"Pixelav Py (Qbin>0, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Prob_y", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[15] = new TH1F("h112",title,100,0.0,1.0);
	sprintf(title,"Pixelav log10(Px) (All Q, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); log10(P_x)", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[17] = new TH1F("h114",title,120,-6.0,0.0);
	sprintf(title,"Pixelav log10(Px) (Qbin>0, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); log10(P_x)", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[19] = new TH1F("h116",title,120,-6.0,0.0);
	sprintf(title,"Pixelav log10(Py) (All Q, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); log10(P_y)", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[21] = new TH1F("h118",title,120,-6.0,0.0);
	sprintf(title,"Pixelav log10(Py) (Qbin>0, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); log10(P_y)", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[23] = new TH1F("h120",title,120,-6.0,0.0);
	sprintf(title,"Pixelav Cluster Charge (size>1, qbin<%d, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Cluster Charge (e)", qbcut, pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[25] = new TH1F("h706",title,nx,0.,300000.);      
	sprintf(title,"Pixelav Cluster Charge (size=1, qbin<%d, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Cluster Charge (e)", qbcut, pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[27] = new TH1F("h708",title,75,0.,75000.);
	sprintf(title,"Normalized Pixelav Cluster Charge (size>1, qbin<%d, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Normalized Cluster Charge (e)", qbcut, pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[29] = new TH1F("h710",title,nx,0.,80000.);      
	sprintf(title,"Pixelav PQ (size=1, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Prob_Q", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[31] = new TH1F("h122",title,100,0.0,1.0);
	sprintf(title,"Pixelav PQ (size>1, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Prob_Q", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[33] = new TH1F("h124",title,100,0.0,1.0);
	sprintf(title,"Pixelav log10(PQ) (size=1, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); log10(P_Q)", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[35] = new TH1F("h126",title,100,-5.0,0.0);
	sprintf(title,"Pixelav log10(PQ) (size>1, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); log10(P_Q)", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[37] = new TH1F("h128",title,100,-5.0,0.0);
	sprintf(title,"Pixelav sizex (p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); sizex (pix)", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[39] = new TH1F("h130",title,10,0.5,10.5);
	sprintf(title,"Pixelav sizey (p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); sizey (pix)", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[41] = new TH1F("h132",title,20,0.5,20.5);
	sprintf(title,"Normalized Pixelav Cluster Charge (all sizes, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Normalized Cluster Charge (e)", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[43] = new TH1F("h712",title,nx,0.,80000.);      
	sprintf(title,"Pixelav sizex (qbin>0, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); sizex (pix)", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[45] = new TH1F("h134",title,10,0.5,10.5);
	sprintf(title,"Pixelav sizey (qbin>0, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); sizey (pix)", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[47] = new TH1F("h136",title,20,0.5,20.5);
 	sprintf(title,"Pixelav Pixel Charge (all clust, p*100>%6.4lf, qscale=%6.3lf, noise=%5.0fe, thresh=%5.0fe, fclus=%5.3f, fgain=%5.3f, rnoise=%5.0fe, lin = %1.1d); Pixel Charge (e)", pcut*100., qscale, noise, q100, common_frac, gain_frac, readout_noise, non_linear);
	hp[48] = new TH1F("h820",title,nx,0.,30000.);
  
   
    printf("mc file number = %d \n", nfile);
	
	// Set style for the the histograms	
	
	hp[1]->SetFillColor(kGreen);
	hp[3]->SetFillColor(kGreen);
	hp[5]->SetFillColor(kGreen);
	hp[7]->SetFillColor(kGreen);
	hp[9]->SetFillColor(kGreen);
	hp[11]->SetFillColor(kGreen);
	hp[13]->SetFillColor(kGreen);
	hp[15]->SetFillColor(kGreen);
	hp[17]->SetFillColor(kGreen);
	hp[19]->SetFillColor(kGreen);
	hp[21]->SetFillColor(kGreen);
	hp[23]->SetFillColor(kGreen);
	hp[25]->SetFillColor(kGreen);
	hp[27]->SetFillColor(kGreen);
	hp[29]->SetFillColor(kGreen);
	hp[31]->SetFillColor(kGreen);
	hp[33]->SetFillColor(kGreen);
	hp[35]->SetFillColor(kGreen);
	hp[37]->SetFillColor(kGreen);
	hp[39]->SetFillColor(kGreen);
	hp[41]->SetFillColor(kGreen);
	hp[43]->SetFillColor(kGreen);
	hp[45]->SetFillColor(kGreen);
	hp[47]->SetFillColor(kGreen);
   hp[48]->SetFillColor(kGreen);

	
	speed = -2;

	iyd = -1; ixd = -1;
	printf("y/x double = %d/%d\n", iyd,ixd);
	
	nevent = 0; ngood = 0; ngood2 = 0; nbad = 0;
	
	ndcol = TYSIZE/2 + 1;
	std::vector<int> ndhit(ndcol, 0);
	
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
		 
		tkp = (float)sqrt((double)(pvec[3]*pvec[3] + pvec[4]*pvec[4] + pvec[5]*pvec[5]));
		if(tkp < tkpmin) continue;
	   ++nevent;
	   
		 // Add noise and analog response to cluster, reformat for flipped barrel coordinate system 
		 
		 triplg(vgauss);
		 for(i=0; i<ndcol; ++i) {ndhit[i] = 0;}
		 icol = 0;
		 if(vgauss[1] < 0.) {icol = 1;}
       pixlst.resize(0);
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
				 }
			 }
		 }
		 
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
        
        
        // Simulate the seed finding
        
        qmax = 0.;
        for(j=0; j<TXSIZE; ++j) {
            for(i=0; i<TYSIZE; ++i) {
                if(clust[j][i] > qmax) {
                    qmax = clust[j][i];
                    max.first = j; max.second = i;
                    
                }
            }
        }
        
		if(qmax < 4000.) continue;
		 		 
	   
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
		qtotal = 0.;
	   int npix = 0;
	   imin = TYSIZE; imax = -1;
	   jmin = TXSIZE; jmax = -1;
       pixIter = pixlst.begin();
	   pixEnd = pixlst.end();
       for ( ; pixIter != pixEnd; ++pixIter ) {
	      j = pixIter->first; 
		  i = pixIter->second;
		  cluster[j][i] = clust[j][i];
		  qtotal += cluster[j][i];
        hp[48]->Fill((double)cluster[j][i]);
		  ++npix;
		   if(j < jmin) {jmin = j;}
		   if(j > jmax) {jmax = j;}
		   if(i < imin) {imin = i;}
		   if(i > imax) {imax = i;}
	   }
		sizex = jmax-jmin+1;
		sizey = imax-imin+1;

// Calculate the hit coordinates in the flipped coordinate system 

	   yhit = -(pvec[0] + (zcen-pvec[2])*pvec[3]/pvec[5]);
	   xhit = -(pvec[1] + (zcen-pvec[2])*pvec[4]/pvec[5]);
	   
// Do the template analysis on the cluster 

       cotalpha = pvec[4]/pvec[5];
	   cotbeta = pvec[3]/pvec[5];
	   eta = fabs(-log((double)(-cotbeta+sqrt((double)(1.+cotbeta*cotbeta)))));
	   
// Combine two single pixels into a double pixel 

	   for(i=0; i<TYSIZE; ++i) {
	      ydouble[i] = false;
	   }
	   if(iyd >= 0 && iyd < TYSIZE) {
	      ydouble[iyd] = true;
		  for(j=0; j<TXSIZE; ++j) {
		     sigtmp = cluster[j][iyd]+cluster[j][iyd+1];
			 cluster[j][iyd] = sigtmp;
		  }
		  for(i=iyd+1; i<TYSIZE-1; ++i) {
		     for(j=0; j<TXSIZE; ++j) {
		        cluster[j][i] = cluster[j][i+1];
		     }
		  }    
		  for(j=0; j<TXSIZE; ++j) {
			 cluster[j][TYSIZE-1] = 0.;
		  }
	   }
		  
	   for(j=0; j<TXSIZE; ++j) {
	      xdouble[j] = false;
	   }
	   if(ixd >= 0 && ixd < TXSIZE-1) {
	      xdouble[ixd] = true;
		  for(i=0; i<TYSIZE; ++i) {
		     sigtmp = cluster[ixd][i]+cluster[ixd+1][i];
			 cluster[ixd][i] = sigtmp;
		  }
		  for(j=ixd+1; j<TXSIZE-1; ++j) {
		     for(i=0; i<TYSIZE; ++i) {
		        cluster[j][i] = cluster[j+1][i];
		     }
		  }    
		  for(i=0; i<TYSIZE; ++i) {
			 cluster[TXSIZE-1][i] = 0.;
		  }
	   }
		
	   
// Do the template analysis on the cluster 
       SiPixelTemplateReco::ClusMatrix clusterPayload{&cluster[0][0], xdouble, ydouble, mrow,mcol};
	   locBx = 1.;
       if(cotbeta < 0.) locBx = -1.;
       locBz = locBx;
       if(cotalpha < 0.) locBz = -locBx;
       ierr = PixelTempReco1D(ID, cotalpha, cotbeta, locBz, locBx, clusterPayload, templ, yrec, sigmay, proby, xrec, sigmax, probx, qbin, speed, probQH);		 

		if(ierr != 0) {
			++nbad;
			if(nbad < 10) {printf("reconstruction failed with error %d \n", ierr);}
		} else {
			++nd1[2];
			if(probQH > 0.5) {
				probQ = 2.*(1. - probQH);
			} else {
				probQ = 2.*probQH;
			}
			if(sizex > 1 || sizey > 1) {++nd1[3];}
			if(probx < 1.e-6) probx = 1.e-6;
			if(proby < 1.e-6) proby = 1.e-6;
			proba = (double)(probx*proby);
			probb = proba*(1.-log(proba));
			log10px=log10((double)probx);
			log10py=log10((double)proby);
			log10pQ=log10((double)probQ);
			qnorm = qtotal/sqrt((double)(1.+cotbeta*cotbeta+cotalpha*cotalpha)); 
			hp[43]->Fill(qnorm);
			if(qbin < qbcut) {
				++nd2[2];
				if(sizex > 1 || sizey > 1) {++nd2[3];}
			}
			if(probb > pcut) {
				++ngood; 
				hp[1]->Fill(qtotal);
				hp[5]->Fill((double)qbin);
				hp[9]->Fill((double)probx);
				hp[13]->Fill((double)proby);
				hp[31]->Fill((double)probQ);
				hp[17]->Fill(log10px);
				hp[21]->Fill(log10py);
				hp[35]->Fill(log10pQ);
				if(sizex == 1 && sizey == 1) {
				   if(qbin<qbcut) {hp[27]->Fill(qtotal);}
				}
				++nd3[2];
				hp[39]->Fill((double)sizex);
				hp[41]->Fill((double)sizey);
				if(qbin > 0) {
				  hp[45]->Fill((double)sizex);
				  hp[47]->Fill((double)sizey);
				}
				if(qbin > 0) {
					hp[11]->Fill((double)probx);
					hp[15]->Fill((double)proby);
					hp[33]->Fill((double)probQ);
					hp[19]->Fill(log10px);
					hp[23]->Fill(log10py);
					hp[37]->Fill(log10pQ);
				}
				if(sizex > 1 || sizey > 1) {
				   ++nd3[3];
				   hp[7]->Fill((double)qbin);			
					qnorm = qtotal/sqrt((double)(1.+cotbeta*cotbeta+cotalpha*cotalpha)); 
					if(qbin < qbcut) {hp[29]->Fill(qnorm);}
				    if(qbin < qbcut && qbin > 0) {hp[3]->Fill(qtotal); ++ngood2;}
					if(qbin<qbcut) {hp[25]->Fill(qtotal);}
				}
				if(qbin < qbcut) {
					++nd4[2];
					if(sizex > 1 || sizey > 1) {++nd4[3];}
				}
			}
			
		}
		if(ngood2 == ntrack2) break;
			
   }
   
/*  Determine current time */

   gettimeofday(&now1, &timz);
   deltas = now1.tv_sec - now0.tv_sec;
   deltaus = now1.tv_usec - now0.tv_usec;
   deltat = ((double)deltaus)/1000000.;
   deltat += (double)deltas;
   printf("ellapsed time = %f seconds \n", deltat);
   
   printf(" total events = %d, template successes = %d, n>1 successes = %d, template failures = %d \n", nevent, ngood, ngood2, nbad);
   printf(" total s=1, s>1 data, s=1, s>1 pixelav = %d  %d  %d  %d\n", nd1[0]-nd1[1], nd1[1], nd1[2]-nd1[3], nd1[3]);
   printf(" qbin<%d s=1, s>1 data, s=1, s>1 pixelav = %d  %d  %d  %d\n", qbcut, nd2[0]-nd2[1], nd2[1], nd2[2]-nd2[3], nd2[3]);
   printf(" prob cut s=1, s>1 data, s=1, s>1 pixelav = %d  %d  %d  %d\n", nd3[0]-nd3[1], nd3[1], nd3[2]-nd3[3], nd3[3]);
   printf(" prob cut, qbin<%d s=1, s>1 data, s=1, s>1 pixelav = %d  %d  %d  %d\n", qbcut, nd4[0]-nd4[1], nd4[1], nd4[2]-nd4[3], nd4[3]);
	
// Create an output filename for this run 
	
//	for(i=0; i<4; ++i) {hp[i]->Fit("landau");}	
	
	sprintf(outfile0,"data_mc_plots/q_dist_pt_bpix%5.5d.pdf[",nfile);
	sprintf(outfile1,"data_mc_plots/q_dist_pt_bpix%5.5d.pdf",nfile);
	sprintf(outfile2,"data_mc_plots/q_dist_pt_bpix%5.5d.pdf]",nfile);
	TCanvas c1("c1", "Charge distributions");
	c1.SetFillStyle(4000);
	c1.Print(outfile0);
	for(i=0; i<16; ++i) {
		hp[i]->Draw();
		c1.Print(outfile1);
	}
   hp[48]->Draw();
   c1.SaveAs("data_mc_plots/pixel_charge.C");
   c1.Print(outfile1);
	hp[3]->Draw();
	hp[2]->Draw("same");
	c1.Print(outfile1);
	gPad->SetLogy(1);
	hp[3]->Draw();
	hp[2]->Draw("same");
	c1.Print(outfile1);
	gPad->SetLogy(1);
	hp[25]->Draw();
	hp[24]->Draw("same");
	c1.Print(outfile1);
	gPad->SetLogy(1);
	hp[27]->Draw();
	hp[26]->Draw("same");
	c1.Print(outfile1);
	gPad->SetLogy(0);
	hp[29]->Draw();
	hp[28]->Draw("same");
	c1.Print(outfile1);
	hp[43]->Draw();
	hp[42]->Draw("same");
	c1.Print(outfile1);
	c1.Clear();
	c1.Divide(1,2);
	c1.Update();
	c1.cd(1);
	hp[1]->Draw();
	hp[0]->Draw("same");
	c1.cd(2);
	gPad->SetLogy(1);
	hp[1]->Draw();
	hp[0]->Draw("same");
	c1.Print(outfile1);
	c1.cd(1);
	gPad->SetLogy(0);
	hp[25]->Draw();
	hp[24]->Draw("same");
	c1.cd(2);
	gPad->SetLogy(1);
	hp[25]->Draw();
	hp[24]->Draw("same");
	c1.Print(outfile1);
	c1.Clear();
	c1.Divide(2,1);
	c1.Update();
	c1.cd(1);
	hp[5]->Draw();
	hp[4]->Draw("same");
	c1.cd(2);
	gPad->SetLogy(1);
	hp[5]->Draw();
	hp[4]->Draw("same");
	c1.Print(outfile1);
	c1.cd(1);
	gPad->SetLogy(0);
	hp[7]->Draw();
	hp[6]->Draw("same");
	c1.cd(2);
	gPad->SetLogy(1);
	hp[7]->Draw();
	hp[6]->Draw("same");
	c1.Print(outfile1);
	c1.Clear();
	c1.Divide(2,1);
	c1.Update();
	c1.cd(1);
	gPad->SetLogy(1);
	hp[9]->Draw();
	hp[8]->Draw("same");
	c1.cd(2);
	gPad->SetLogy(1);
	hp[13]->Draw();
	hp[12]->Draw("same");
	c1.Print(outfile1);
    c1.cd(1);
    gPad->SetLogy(1);
    hp[17]->Draw();
    hp[16]->Draw("same");
    c1.cd(2);
    gPad->SetLogy(1);
    hp[21]->Draw();
    hp[20]->Draw("same");
    c1.Print(outfile1);
	c1.cd(1);
	gPad->SetLogy(1);
	hp[11]->Draw();
	hp[10]->Draw("same");
	c1.cd(2);
	gPad->SetLogy(1);
	hp[15]->Draw();
	hp[14]->Draw("same");
	c1.Print(outfile1);
    c1.cd(1);
    gPad->SetLogy(1);
    hp[19]->Draw();
    hp[18]->Draw("same");
    c1.cd(2);
    gPad->SetLogy(1);
    hp[23]->Draw();
    hp[22]->Draw("same");
    c1.Print(outfile1);
	c1.cd(1);
	gPad->SetLogy(1);
	hp[31]->Draw();
	hp[30]->Draw("same");
	c1.cd(2);
	gPad->SetLogy(1);
	hp[33]->Draw();
	hp[32]->Draw("same");
	c1.Print(outfile1);
	c1.cd(1);
	gPad->SetLogy(1);
	hp[35]->Draw();
	hp[34]->Draw("same");
	c1.cd(2);
	gPad->SetLogy(1);
	hp[37]->Draw();
	hp[36]->Draw("same");
	c1.Print(outfile1);
	c1.cd(1);
	gPad->SetLogy(1);
	hp[31]->Draw();
	hp[30]->Draw("same");
	c1.cd(2);
	gPad->SetLogy(1);
	hp[35]->Draw();
	hp[34]->Draw("same");
	c1.Print(outfile1);
	c1.cd(1);
	gPad->SetLogy(0);
	hp[39]->Draw();
	hp[38]->Draw("same");
	c1.cd(2);
	gPad->SetLogy(0);
	hp[41]->Draw();
	hp[40]->Draw("same");
	c1.Print(outfile1);
	c1.cd(1);
	gPad->SetLogy(0);
	hp[45]->Draw();
	hp[44]->Draw("same");
	c1.cd(2);
	gPad->SetLogy(0);
	hp[47]->Draw();
	hp[46]->Draw("same");
	c1.Print(outfile1);
	
	//   THStack hs0("hs0","Cluster Charge (electrons)");
	//   hs0.Add(hp[0]);
	//   hs0.Add(hp[1]);
	//   hs0.Draw();
	//   c1.Print(outfile1);
	//   THStack hs1("hs1","Normalized Cluster Charge (electrons)");
	//   hs1.Add(hp[2]);
	//   hs1.Add(hp[3]);
	//   hs1.Draw();
	//   c1.Print(outfile1);
	c1.Print(outfile2);
	//   TFile f("pixelav_clustx_vs_cotalpha.root","new");
	//   pp[3]->Write();
	//   f.ls();
	//   f.Close();
	
return 0;
} // MAIN__ 




// ***************************************************************** 
//! Calculate 21 gaussianly-distributed random numbers.
//! \param x - a vector holding 21 random numbers generated with rms = 1.            
// ***************************************************************** 
int triplg(std::vector<float>& x)
{
    // Initialized data 

    static int fclusl = -1;

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

    if (fclusl) {
	twopi = 2.*acos((double)-1.);
	ibase = TYTEN;
	fclusl = 0;
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



// ***************************************************************** 
//! Calculate 21 uniformly-distributed random numbers from -1 to +1.
//! \param x - a vector holding 21 random numbers.            
// ***************************************************************** 
int triplu(std::vector<float>& x)
{
    // Initialized data 
	
    // Local variables 
    static int i__;
	
    // Function Body 
	
	//  Initalize the parameters 
	
    for (i__ = 0; i__ < TYSIZE; ++i__) {
		x[i__] = 2.*((float)random())/((float)RAND_MAX)-1.;
    }
    return 0;
} // triplu 

