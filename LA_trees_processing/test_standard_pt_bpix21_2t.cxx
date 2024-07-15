//! \file template_code6d.cpp
//!
//! Apply standard CMSSW pixel reconstruction to pixelav hits
//! Do Lorentz angle extraction using cluster summary lorentz_clust.txt and the pixel input file.
//! Update to new summary file format
//! Add second threshold
//! All layer info
//! Take simulation charge scale from template
//! Add Plus/Minus end of layers
//! Update for Phase 1 vcal values
//! Add flipped vs unflipped modules

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

Double_t fitf(Double_t *x,Double_t *par)
{
Double_t arg;
  if(x[0] < par[3]) {
     arg = par[1]*par[1] + par[2]*par[2]*(x[0]-par[3])*(x[0]-par[3]); 
  } else {
     arg = par[1]*par[1] + par[4]*par[4]*(x[0]-par[3])*(x[0]-par[3]); 
  }  
  Double_t fitval = par[0]+sqrt(arg); 
  return fitval;
  
}
Double_t singlef(Double_t *x,Double_t *par)
{
	Double_t arg = par[1]*par[1] + par[2]*par[2]*(x[0]-par[3])*(x[0]-par[3]); 

	Double_t fitval = par[0]+sqrt(arg); 
	return fitval;
	
}

// Main program  

int main(int argc, char *argv[])
{
   
    std::vector<float> pvec(6), wgauss(TYSIZE), vgauss(TYSIZE), xgauss(TYSIZE), ygauss(TYSIZE), zgauss(TYSIZE);
    static bool fpix;
    float pixin[TXSIZE][TYSIZE];
    bool ydouble[TYSIZE], xdouble[TXSIZE];
    static float thick, xsize, ysize, noise, zcen, xcmssw, ycmssw, gain_frac, readout_noise, q100_frac, common_frac;
	static float xhit, yhit, sigmax, sigmay, probx, proby, signal, cotalpha, cotbeta, probQH, probQ, tkp, tkpmin, probmin, probxy, probQmin, chimax, qscale, locBx, locBz;
    static int sizex, sizey, icol, ndcol;
    static int ndata, nfile, neh, nevent, ID0, ID, non_linear, linecolor, layer, end, ltpt; 
	static vector<int> nbin(5,0);
	int mrow = TXSIZE, mcol = TYSIZE;
    int i, j, k, ierr, qbin, ntrack, jmin, jmax, imin, imax, numadd, flipped, module, idcol, lay, flp;
	double dx, dy, eta, dxc, dyc, qtotal, adc;
	const double gain = 3.19;
	const double ped = 16.46;
	const double p0 = 0.01218;
	const double p1 = 0.711;
	const double p2 = 203.;
	const double p3 = 148.;	
	static double vcal = 47.;	
	static double vcaloffst = 60.;
    static float qseed = 5000.;
	static int iyd, ixd, speed;	
	static float q100, q101, q50, q10, qmax; 
	float sigtmp, qin, chisq;
    static char infile[80], header[80], c, outfile0[80], outfile1[80], outfile2[80];
	float rocsim(float, int, bool);
	int triplg(std::vector<float>&);
//	int random(void);
    Double_t fitf(Double_t *,Double_t *);
	float cluster[TXSIZE][TYSIZE], clust[TXSIZE][TYSIZE], rclust[TXSIZE][TYSIZE];
	bool bclust[TXSIZE][TYSIZE];
    std::pair<int, int> pixel, max;

    FILE *ifp;
	
   struct timeval now0, now1;
   struct timezone timz;
   long deltas, deltaus;
   double deltat;
	
	// Ask for external input into which pixels to join 
	
	printf("enter momentum cut (GeV), probxy_min, probQ_min, ID0, layer (1-6,0=all,7=all new), Min/Pls/Both (0/1/2), unf/flp/both (0/1/2), lt/pt (0/1) \n");
	scanf("%f %f %f %d %d %d %d %d", &tkpmin, &probmin, &probQmin, &ID0, &layer, &end, &flp, &ltpt);
	printf("momentum cut = %f, prob_min = %f, probQ_min = %f, ID0 = %d, layer = %d, end = %d, flipped = %d, lt/pt = %d \n", tkpmin, probmin, probQmin, ID0, layer, end, flp, ltpt);
	linecolor = 2;
    
	double  halfxs=300.;
	int nx=120;	
	gStyle->SetOptStat(1111);
	gStyle->SetOptFit(1111);
	gStyle->SetHistLineWidth(2);
	static vector<TH1F*> hp(9);
	hp[0] = new TH1F("h101","dy (all sig); #Deltay (#mum)",nx,-halfxs,halfxs);
	hp[1] = new TH1F("h102","dx (all sig); #Deltax (#mum)",nx,-halfxs,halfxs);
	hp[2] = new TH1F("h300","cotbeta (probx<10-3)",nx,-10.,10.);      
	hp[3] = new TH1F("h301","cotalpha (probx<10-3)",nx,-0.25,0.25);   
	hp[4] = new TH1F("h103","dx (size = 1); #Deltax (#mum)",nx,-halfxs,halfxs);
	hp[5] = new TH1F("h104","dx_temp (size > 1); #Deltax (#mum)",nx,-halfxs,halfxs);
	hp[6] = new TH1F("h105","data n_pix = 1 clusters; cot(#alpha)",50,-1.5,1.0);
	hp[7] = new TH1F("h106","sim n_pix = 1 clusters; cot(#alpha)",50,-1.5,1.0);
	hp[8] = new TH1F("h107","chisquare; #chi^2",50,0.,25.);

	// Set style for the the histograms	
	
	for(i=0; i<7; ++i) {
		hp[i]->SetLineColor(2);
		hp[i]->SetFillColor(38);
	}
	hp[7]->SetLineColor(4);
	hp[8]->SetLineColor(2);
	hp[8]->SetFillColor(38);
	
	// Make some profile Histograms
	
	static vector<TProfile*> pp(10);
	pp[0] = new TProfile("pqhycms","y_cmssw",11,0,2.75,"s"); 
	pp[1] = new TProfile("pqhxcms","x_cmssw",11,0,2.75,"s"); 
	pp[2] = new TProfile("pixnycotb","ysize(pix); cot(#beta)",40,-4.,4.," "); 
	pp[3] = new TProfile("pixnxcota","xsize(pix); cot(#alpha)",40,-0.80,0.80," "); 
	pp[4] = new TProfile("datnxcota","xsize(pix); cot(#alpha)",40,-0.80,0.80," "); 
	pp[5] = new TProfile("datnycotb","ysize(pix); cot(#beta)",40,-4.,4.," "); 
	pp[6] = new TProfile("pixnxcotany1","xsize(pix) [ny=1]; cot(#alpha)",10,-0.20,0.20," "); 
	pp[7] = new TProfile("datnxcotany1","xsize(pix) [ny=1]; cot(#alpha)",10,-0.20,0.20," "); 
	pp[8] = new TProfile("datnycotb14","ysize(pix); cot(#beta)",40,-4.,4.," "); 
	pp[9] = new TProfile("datnycotb58","ysize(pix); cot(#beta)",40,-4.,4.," "); 
	
	// Set style for the the profile histograms	
	
	for(i=0; i<10; ++i) {
		pp[i]->SetLineColor(4);
//		pp[i]->SetStats(kFALSE);
		pp[i]->SetMinimum(0.);
	}
	pp[2]->SetLineColor(linecolor);
	pp[2]->SetLineStyle(2);
	pp[3]->SetLineColor(linecolor);
	pp[3]->SetLineStyle(2);
	pp[6]->SetLineColor(linecolor);
	pp[6]->SetLineStyle(2);
	
	
//  Read which data and inputs to use (use c file i/o which is much faster than c++ i/o) 
	
	ifp = fopen("q_dist_2t.txt", "r");
	if (ifp==NULL) {
      printf("no pixel initialization file found/n");
      return 0;
	}
	
	fscanf(ifp,"%d %d %f %f %f %f %f %f %f %d", &ndata, &nfile, &noise, &q100, &q101, &q100_frac, &common_frac, &gain_frac, &readout_noise, &non_linear);
	fclose(ifp);
	printf("data file %d, mc file %d noise = %f, threshold = %f, threshold1 = %f, rms threshold frac = %f, common_frac = %f, gain fraction = %f, readout noise = %f, nonlinear_resp = %d \n", ndata, nfile, noise, q100, q101, q100_frac, common_frac, gain_frac, readout_noise, non_linear);
		
//  Create an input data file name for this run 
	
	if(ltpt == 0) {sprintf(infile,"Efield_output/lt_bpix%6.6d.txt",ndata);} else {sprintf(infile,"pt_bpix%6.6d.txt",ndata);}
	
//  Open input file and read header info 
	
    ifp = fopen(infile, "r");
    if (ifp==NULL) {
		printf("no pixel data file found/n");
		return 0;
    }
	
	ntrack = 0;
    while(fscanf(ifp,"%f %f %lf %d %d %d %f %f %lf %lf %f %f %d %d %d", &cotalpha, &cotbeta, &qtotal, &sizex, &sizey, &qbin, &probx, &proby, &dx, &dy, &probQH, &tkp, &flipped, &module, &lay) != EOF) {
		 if(probQH > 0.5) {
			 probQ = 2.*(1. - probQH);
		 } else {
			 probQ = 2.*probQH;
		 }
//      if(layer == 7 && lay > 3) lay = 7;
		if(tkp < tkpmin) continue;
        if(layer > 0 && lay != layer) continue;
        if(end == 0 && module > 0) continue;
        if(end == 1 && module < 0) continue;
        if(flp < 2 && flp != flipped) continue;
		probxy = probx*proby*(1. - (float)log((double)(probx*proby)));
		if(probxy < probmin) continue;
		if(probQ < probQmin) continue;
		++ntrack; 
//		if(sizex ==1 && sizey == 1) continue;
		if(qbin > 0 && qbin < 4) {pp[4]->Fill((double)cotalpha, (double)sizex);}
		if(qbin > 0 && qbin < 4) {
            pp[5]->Fill((double)cotbeta, (double)sizey);
            if(module < 0) {pp[8]->Fill((double)cotbeta, (double)sizey);} else {pp[9]->Fill((double)cotbeta, (double)sizey);}            
        }
		if(qbin > 0 && qbin < 4 && sizey == 1) {pp[7]->Fill((double)cotalpha, (double)sizex);}
	}
    fclose(ifp);
	
	printf("number of tracks = %d \n", ntrack);

//  Create an input filename for this run 

    if(nfile < 10000) {
        
        sprintf(infile,"simulated_clusters/template_events_d%4.4d.out",nfile);
        
    } else {
        
        sprintf(infile,"simulated_clusters/template_events_d%5.5d.out",nfile);
        
    }

//  Open input file and read header info 

	ifp = fopen(infile, "r");
    if (ifp==NULL) {
		printf("no pixel simulated data file '%s' found\n", infile);
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
	
	
	static vector<float> sxc(4,0.), sxc2(4,0.), syc(4,0.), syc2(4,0.); 

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
	
	iyd = -1; ixd = -1;
	printf("y/x double = %d/%d\n", iyd,ixd);
	
	int totale=0;  int goodx=0; int goody=0;
	
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
		 
       cotalpha = pvec[4]/pvec[5];
		 cotbeta = pvec[3]/pvec[5];
		 eta = fabs(-log((double)(-cotbeta+sqrt((double)(1.+cotbeta*cotbeta)))));
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
        
		if(qmax < qseed) continue;
		 
	   
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
		imin = TYSIZE; imax = 0;
		jmin = TXSIZE; jmax = 0;
       for ( ; pixIter != pixEnd; ++pixIter ) {
	      j = pixIter->first; 
		  i = pixIter->second;
		  cluster[j][i] = clust[j][i];
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

		 speed = -2;
	   
// Do the template analysis on the cluster 
       SiPixelTemplateReco::ClusMatrix clusterPayload{&cluster[0][0], xdouble, ydouble, mrow,mcol};
	   locBx = 1.;
       if(cotbeta < 0.) locBx = -1.;
       locBz = locBx;
       if(cotalpha < 0.) locBz = -locBx;
       ierr = PixelTempReco1D(ID, cotalpha, cotbeta, locBz, locBx, clusterPayload, templ, ycmssw, sigmay, proby, xcmssw, sigmax, probx, qbin, speed, probQ);		 

	   if(ierr != 0) {
	      printf("reconstruction failed with error %d \n", ierr);
		} else {
	   
	      ++nbin[0];
// Check resolution and weights 
		   if(iyd != 0) {
		     dyc = ycmssw - (TYSIZE/2)*ysize - yhit;
		   } else {
		     dyc = ycmssw - ((TYSIZE/2)-0.5)*ysize - yhit;
		   }
		   syc[0] += dyc; syc2[0] += dyc*dyc; 
		   if(ixd != 0) {
			 dxc = xcmssw - (TXSIZE/2)*xsize - xhit;
		   } else {
			 dxc = xcmssw - ((TXSIZE/2)-0.5)*xsize - xhit;
		   }
			chisq = dyc*dyc/(sigmay*sigmay) + dxc*dxc/(sigmax*sigmax);
			hp[8]->Fill((double)chisq);
			if(probQH > 0.5) {
				probQ = 2.*(1. - probQH);
			} else {
				probQ = 2.*probQH;
			}
			probxy = probx*proby*(1. - (float)log((double)(probx*proby)));
			if(probxy < probmin) continue;
			if(probQ < probQmin) continue;
		   sxc[0] += dxc; sxc2[0] += dxc*dxc; 
		   hp[2]->Fill((double)cotbeta);
		   hp[3]->Fill((double)cotalpha);
		   if(sizex ==1 && sizey == 1) {hp[7]->Fill((double)cotalpha);}
		   
		   hp[0]->Fill(dyc);
		   hp[1]->Fill(dxc);
		   if(sizex == 1) {hp[4]->Fill(dxc);} else {hp[5]->Fill(dxc);}
//			if(chisq > chimax) continue;
		   pp[0]->Fill(eta, dyc);
		   pp[1]->Fill(eta, dxc);
//		   pp[2]->Fill((double)cotbeta,(double)sizey);
//			if(sizex == 1 && sizey == 1) continue;
		   if(qbin > 0 && qbin < 4) {pp[3]->Fill((double)cotalpha,(double)sizex);}
		   if(qbin > 0 && qbin < 4) {pp[2]->Fill((double)cotbeta, (double)sizey);}
			if(qbin > 0 && qbin < 4 && sizey == 1) {pp[6]->Fill((double)cotalpha, (double)sizex);}
//		   h2[0]->Fill((double)cotbeta,(double)cotalpha);

	    }
		
   }
   
/*  Determine current time */

   gettimeofday(&now1, &timz);
   deltas = now1.tv_sec - now0.tv_sec;
   deltaus = now1.tv_usec - now0.tv_usec;
   deltat = ((double)deltaus)/1000000.;
   deltat += (double)deltas;
   printf("ellapsed time = %f seconds \n", deltat);
   
   printf(" total events = %d, probx > 10^{-3} = %d, proby > 10^{-3} = %d \n", totale, goodx, goody);
   printf(" low q failures = %d \n", nbin[4]);
	   
   printf(" CMSSW algorithm \n");
	   
   for(j=0; j<1; ++j) {
      syc[j] /= (float)nbin[j]; syc2[j] /= (float)nbin[j];
      syc2[j] = sqrt((double)(syc2[j] - syc[j]*syc[j]));
      printf(" avg y residual[%1d] = %f +- %f \n", j, syc[j], syc2[j]);       
      sxc[j] /= (float)nbin[j]; sxc2[j] /= (float)nbin[j];
      sxc2[j] = sqrt((double)(sxc2[j] - sxc[j]*sxc[j]));
      printf(" avg x residual[%1d] = %f +- %f \n", j, sxc[j], sxc2[j]); 
  }
  
// Make resolution graphs from profile errors

  float binx, biny;
  for(i=0; i<2; ++i) {
  printf("\n Profile errors %d, %s \n",i, pp[i]->GetTitle());
     for(j=0; j<11; ++j) {
	    binx=pp[i]->GetBinCenter(j+1); biny=pp[i]->GetBinContent(j+1);
        printf("%f  %f \n",binx, biny);
	 }
  }
  for(i=2; i<4; ++i) {
  printf("\n Profile errors %d, %s \n",i, pp[i]->GetTitle());
     for(j=0; j<80; ++j) {
	    binx=pp[i]->GetBinCenter(j+1); biny=pp[i]->GetBinContent(j+1);
        printf("%f  %f \n",binx, biny);
	 }
  }
  
   TCanvas c1("c1", header);
   c1.SetFillStyle(4000);
/*
 * Histograms plotting
 */
   for(i=0; i<2; ++i) {hp[i]->Fit("gaus"); hp[i+4]->Fit("gaus");}
   
   TF1 *func = new TF1("func", singlef, -0.50, 0.40, 4);
   func->SetParameters(1.,0.1,1.6,-0.4,1.2);
   func->SetParNames ("Offset","RMS Constant","SlopeL","cot(alpha)_min");
	func->SetLineColor(4);
	func->SetLineStyle(2);
	func->SetLineWidth(1);
	pp[3]->Fit(func, "R");
//	pp[3]->Draw();
//	c1.SaveAs("pixelav_clustx_vs_cotalpha.C");
	func->SetLineColor(2);
	pp[4]->Fit(func, "R");
	
// Create an output filename for this run 
	
	sprintf(outfile0,"data_mc_plots/pixel_histos%4.4d.pdf[",nfile);
	sprintf(outfile1,"data_mc_plots/pixel_histos%4.4d.pdf",nfile);
	sprintf(outfile2,"data_mc_plots/pixel_histos%4.4d.pdf]",nfile);
	c1.Print(outfile0);
	for(i=0; i<9; ++i) {
       hp[i]->Draw();
	   c1.Print(outfile1);
	}
	for(i=0; i<6; ++i) {
	   pp[i]->Draw();
	   c1.Print(outfile1);
	}
	pp[3]->SetStats(kFALSE);
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(kFALSE);
	pp[3]->Draw();
	pp[4]->SetStats(kFALSE);
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(kFALSE);
	pp[4]->Draw("same");
	c1.Print(outfile1);
	pp[2]->SetStats(kFALSE);
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(kFALSE);
	pp[2]->Draw();
	pp[5]->SetStats(kFALSE);
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(kFALSE);
	pp[5]->Draw("same");
	c1.Print(outfile1);
	pp[2]->SetStats(kFALSE);
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(kFALSE);
	pp[2]->Draw();
	pp[8]->SetStats(kFALSE);
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(kFALSE);
	pp[8]->Draw("same");
	c1.Print(outfile1);
	pp[2]->SetStats(kFALSE);
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(kFALSE);
	pp[2]->Draw();
	pp[9]->SetStats(kFALSE);
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(kFALSE);
	pp[9]->Draw("same");
	c1.Print(outfile1);
	pp[6]->SetStats(kFALSE);
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(kFALSE);
	pp[6]->Draw();
	pp[7]->SetStats(kFALSE);
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(kFALSE);
	pp[7]->Draw("same");
	c1.Print(outfile1);
	c1.Clear();
	c1.Divide(2,1);
	c1.Update();
	c1.cd(1);
	pp[3]->Draw();
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(kFALSE);
	pp[4]->Draw("same");
	c1.SaveAs("data_mc_plots/pixelav_clustx_vs_cotalpha_bpix.C");
	c1.cd(2);
	pp[2]->Draw();
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(kFALSE);
	pp[5]->Draw("same");
	c1.Print(outfile1);
	
	//   h2[0]->Draw();
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

    static int fcall = -1;

    // Local variables 
    static float r1, r2;
    static int i__;
    static float r__;
    static int ibase;
    static std::vector<float> rbuff(210);
    static float twopi;
    double arg, phi;



    // Function Body 

//  Initalize the parameters 

    if (fcall) {
	twopi = 2.*acos((double)-1.);
	ibase = 210;
	fcall = 0;
    }

//  If all random numbers used up, generate 210 more 

    if (ibase == 210) {
	   for (i__ = 0; i__ < 209; i__ += 2) {
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
    for (i__ = 0; i__ < 21; ++i__) {
	   x[i__] = rbuff[ibase + i__];
    }
    ibase += 21;
    return 0;
} // triplg 





// ***************************************************************** 
//! Simulate the analog electronics for PSI46v2 (10/18/06).
//! Returns digitized information as float.    
//!
//! \param qin   - input signal in electrons               
//! \param iresp - selects from among several response functions     
//!                0 to return CMSSW adc units and linear response
//!                1 PSI46v2 response function R78 (lin at large Q)  
//!                2 PSI46v2 response function CMSSW(lin at med  Q)
//!                3 PSI46v2 response function mod R78 (lin at low Q)
//!                4 PSI46v2 response function R74 (default)
//!                5 PSI46v2 response function R76 (PSI optimized)
//! \param inverse - false: convert charge to adc, true: adc to charge     
// ***************************************************************** 
float rocsim(float qin, int iresp, bool inverse)
{
    // Initialized data 

    const float vcal = 65.;
	const float p0[5] = {0.00373,0.00382,0.00382,0.00388,0.00359};
	const float p1[5] = {1.240,0.886,0.40,1.015,1.067};
	const float p2[5] = {121.5,112.7,120.,104.,124.};
	const float p3[5] = {125.6,113.0,80.,131.,126.};
	double q, adc, arg, tanhin;
	float roc;
	int ih;
	static float roc0;
	static int iresp_cur;
	static bool fcall=true;
    // Function Body 
	
	if(fcall || iresp != iresp_cur) {
	   iresp_cur = iresp;
	   if(iresp < 1) {
	      roc0=0.;
	   } else {
		     ih=iresp-1;
	         roc0=p3[ih]-p2[ih]*tanh((double)p1[ih]);
	   }
	   fcall=false;
	}
	
//  If the signal in adc units is requested

    if(!inverse) {

//  If input charge is negative, set to zero 

       q = (double)qin;
       if (q < 0.) { q = 0.;}
	   if(iresp < 1) {
	      roc = q/135.;
	   } else {
	      ih=iresp-1;
		  if (ih > 4) {ih=4;}
	      roc=p3[ih]+p2[ih]*tanh(p0[ih]*(q/vcal) - p1[ih])-roc0;
	   }
	} else {
	
//  If the inverse signal in charge units is requested

       adc = (double)qin;
       if (adc < 0.) { adc = 0.;}
	   if(iresp < 1) {
	      roc = adc*135.;
	   } else {
	      ih=iresp-1;
		  if (ih > 4) {ih=4;}
		  arg = (adc+roc0-p3[ih])/p2[ih];
		  if (fabs(arg) < 1.0) {
		     tanhin=log((1.+arg)/(1.-arg))/2.;
		     roc=vcal*(tanhin+p1[ih])/p0[ih];
		  } else {
		     roc=40000.;
		  }
	   }
	}
	
	return roc;
} // rocsim 

