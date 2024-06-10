 /* pixel32.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* Fully vectorize the charge deposition and include the effect of the magnetic field on the delta-rays */
/* Use Henrich mobility and Llubjana trapping */
/* Reduce Hall factor to 1.12 (12/10/03) */
/* Adaptive step sizing using Cash-Karp embedded 5th-order technique: version (6/23/05) */
/* Add header output to barrel_ten.out.  Switched on by adding the character '1' to the */
/* beginning of the ascii header in pixel.init (7/05/05) */
/* Add pion energy dependence of cross sections from H. Bichsel.  Use magnitude of pion direction to */
/* store the information.  Assumes 45 GeV if not specified in old files (11/10/05) */
/* Change pixel array to 21x7 to accommodate wider range of input angles (04/10/06) */
/* Version to automatically generate multiple output files while incrementing the cluster length */
/* Add electron hall factor rhe to input list */
/* Add NIST Estar inverse stopping powers: drde (11/15/2007) */
/* Add multiple scattering to primary delta rays for NIST Estar choice (11/14/2009) */
/* Increase event array size to 21x13 (11/21/2008) */
/* Multiple run version for use on cluster */
/* Enforce efield boundary conditions */
/* Fix factor of 2 in diffusion equation */
/* Allow independent adjustment of fluxes to separately tune e/h trapping rates */
/* Needs non-standard pixel2.init initialization files */
/* Include external weighting potential lookup table in wgt_pot.init */
/* Version for 3d and planar sensors */
/* Restructured to use common code for all calculations */
/* Use sse2neon.h to run on arm based processors */

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stddef.h>
#ifdef __POWERPC__
#include <altivec.h>
#undef pixel
#else
#if defined(__arm64__) || defined(__ARM_NEON)
#include "sse2neon.h"
#else
#include <xmmintrin.h>
#endif
#endif

/* Global symbols */

/* Maximum number of e-h pairs to be stored in static arrays */

#define NEHSTORE 500000

/* Define pixel signal buffer sizes */

#define TXSIZE 21
#define TYSIZE 13
/* Define maximum E-field and weighting potential array sizes */

#define NARRAYX 26
#define NARRAYY 26
#define NARRAYZ 94

/* Define the maxumum number of runs */

#define TEMPMAX 500

/* Table of constant values */

static int c__1 = 1;
static int c__2 = 2;
static int c__3 = 3;
static int c__4 = 4;
static int c__120 = 120;

/* Prescaling factor for electrons and holes: transport only 1/Nscale carriers to save time */

static int Nscale = 10;  /* This doesn't cause additional fluctuations (we already get 22,000 e-h pairs per 300um Si) */

   typedef union vect_or_f {
#ifdef __POWERPC__
     vector float v;
#else
     __m128 v;
#endif
     float f[4];
   } vect_or_f;
#ifdef __POWERPC__   
   typedef union vect_or_c {
     vector unsigned char v;
     unsigned char c[16];
   } vect_or_c;
#endif
   typedef union vect_or_i {
#ifdef __POWERPC__   
     vector unsigned int v;
#else
  __m128 v;
#endif
     unsigned int i[4];
   } vect_or_i;

/* Define the E-field array of vectors and associated quantities */

    static vect_or_f efield[NARRAYX][NARRAYY][NARRAYZ];
    static vect_or_f wgtpot[NARRAYX][NARRAYY][NARRAYZ][3];
    static int npixx, npixy, npixz;    
/* Define array to help with bounds checking */
    static int mnode[3];
    static vect_or_f bfield;
    static float bfield_z;
    static char header[80];
    
/* function prototypes */

#include "pixelav3d_prototypes.h"

   int runinit(int *, int *, int *, float *, float *, float *, float *, float *, float *);

/* Main program */ int main(int argc, char *argv[])
{
    /* System generated locals */
    int i__1;
    float r__1, r__2;
    double d__1, d__2;

    /* Local variables */
    static float vect[6];
    static float thick, xsize, ysize, temp, flux[2], rhe, rhh, impwdth, implenn, implenp;    
    static int i__, indeh[2][NEHSTORE]	/* was [2][300000] */;
    static int nto2in, lux, initseed, ivec[25];
    static float pixel[2][TXSIZE][TYSIZE];
    static int ntotin, neh;
    static vect_or_f xeh[2][NEHSTORE]	/* was [4][2][300000] */;
    time_t now;
    struct tm *nows;
    int sec, min, hour, yday, j, k;
    static int fileind, filebase, fileoff, runsize, irun, ievent, frun, nrun, procid, new_drde, ehole, iruny, irunx, nxrun, nyrun;
    static float rvec[4], pimom, xoffset, yoffset, lenxmin, lenxmax, deltaxlen, lenymin, lenymax, deltaylen, locdir[3], cotalpha, cotbeta;
	static float clusxlen, clusylen;
	static char outfile[80], seedfile[80];
	static double alpha;

    FILE *isfp, *iifp, *ofp;

    
	/* If no arguments, quit */
	
    if(argc < 2) {
		printf("Need at least one argument to specify run \n");
		return 0;
    }
    
	/* A single argument is a first run number from the runlist */
	
    if(argc == 2) {
		sscanf(argv[1],"%d", &frun);
		if(frun < 1 || frun > TEMPMAX) {printf("frun %d is illegal, quit \n", frun); return 0;}
		nrun = 1;
		printf("Starting from runlist number %d, processing %d runs\n", frun, nrun);
    }
    
	/* If two arguments, second could be a number of runs or a fork instruction */
	
    if(argc == 3) {
		sscanf(argv[1],"%d", &frun);
		if(frun < 1 || frun > TEMPMAX) {printf("frun %d is illegal, quit \n", frun); return 0;}
		if(*argv[2] == 'f') {nrun = 1;} else {
			sscanf(argv[2],"%d", &nrun);
			if(nrun < 1 || nrun > (TEMPMAX-frun)) {printf("nrun %d is illegal, quit \n", nrun); return 0;}
		}	
		printf("Starting from runlist number %d, processing %d runs\n", frun, nrun);
		if(*argv[2] == 'f') {
			procid = fork();
			if(procid) {
				printf("Forking process, id = %d\n", procid);
				return 0; 
		    }			
		}
	}
	
	/* If three arguments, retrieve first run, number of runs, and possible fork command */
	
	if(argc == 4) {
		sscanf(argv[1],"%d", &frun);
		if(frun < 1 || frun > TEMPMAX) {printf("frun %d is illegal, quit \n", frun); return 0;}
		sscanf(argv[2],"%d", &nrun);
		if(nrun < 1 || nrun > (TEMPMAX-frun)) {printf("nrun %d is illegal, quit \n", nrun); return 0;}
		printf("Starting from runlist number %d, processing %d runs\n", frun, nrun);
		if(*argv[3] == 'f') {
			procid = fork();
			if(procid) {
				printf("Forking process, id = %d\n", procid);
				return 0; 
			}			
	    }
    }
    
	/*  Define the detector parameters from the global initialization files */
	
    if(pixinit(&pimom, &thick, &xsize, &ysize, &impwdth, &implenn, &implenp, &temp, flux, &rhe, &rhh, &ehole, &new_drde, &filebase) == 0) {return 0;}
    
    if(pimom < 1.1) {pimom = 45.;}
	
	/*  Define the run parameters from the local run initialization file */
	
    runinit(&initseed, &fileoff, &runsize, &lenxmin, &lenxmax, &deltaxlen, &lenymin, &lenymax, &deltaylen);
	
	/* the file index is the sum of an overall base number and a local run offset */
	
	fileind = filebase + fileoff + frun;
	
	/*  Create a seedfile name for this run */
	
    sprintf(seedfile,"seedfile%5.5d",fileind);
    
/*  Determine current time */

    now = time(NULL);
    nows = localtime(&now);
    sec = (*nows).tm_sec;
    min = (*nows).tm_min;
    hour = (*nows).tm_hour;
    yday = (*nows).tm_yday;
    
    printf("Begin on day %d at %02d:%02d:%02d\n", yday, hour, min, sec);

/*  Initialize the random number generation */

/* First check to see if any intermediate state has been saved */

    isfp = fopen(seedfile, "r");
    if (isfp==NULL) {

/* If no seedfiled, use single seed and set-up generator */

       lux = 3;
       ntotin = 0;
       nto2in = 0;
       rluxgo_(&lux, &fileind, &ntotin, &nto2in);
	   irun=0; ievent=0;
    } else {
    
/* read-in 25 ints and restore generator to previous state */

       fscanf(isfp,
       "%x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %d %d", 
       &ivec[0], &ivec[1], &ivec[2], &ivec[3], &ivec[4], &ivec[5], &ivec[6], &ivec[7],
       &ivec[8], &ivec[9], &ivec[10], &ivec[11], &ivec[12], &ivec[13], &ivec[14], &ivec[15],
       &ivec[16], &ivec[17], &ivec[18], &ivec[19], &ivec[20], &ivec[21], &ivec[22], &ivec[23],
       &ivec[24], &irun, &ievent);
       fclose(isfp);
       rluxin_(ivec);
    }
    
	
	/* Loop until all runs are complete */
	
    while(irun < nrun) {
		
/*  Create a filename for this run */
	
	   fileind = filebase + fileoff + frun + irun;

       sprintf(outfile,"template_events_d%5.5d.out",fileind);
	
	   if(frun == 1 && irun == 0 && ievent == 0) {
		
		/*  copy  header to the output file on first event */
		
		   ofp = fopen(outfile, "w");
		   fprintf(ofp,"%s \n", &header[0]);
		   fprintf(ofp,"%f  %f  %f\n", xsize, ysize, thick);
		   fclose(ofp); 	   
	   }
	
		   
	   while(ievent < runsize) {
			   
	/* Generate initial position and direction of the track */
			   
		   ranlux_(rvec,&c__4);
		   clusxlen = (lenxmin + rvec[2]*(lenxmax-lenxmin))*xsize;
		   cotbeta = clusxlen/thick;	
		   clusylen = (lenymin + rvec[3]*(lenymax-lenymin))*ysize;
		   cotalpha = clusylen/thick;
		   locdir[2] = 1./sqrt((double)(1.+cotbeta*cotbeta+cotalpha*cotalpha));
		   locdir[0] = cotbeta*locdir[2];
		   locdir[1] = cotalpha*locdir[2];
			   
			   /*  Calculate the offsets from the detector center to its front face */
			   
		   xoffset = locdir[0]/locdir[2] * thick / 2.;
		   yoffset = locdir[1]/locdir[2] * thick / 2.;
			   
		   if(locdir[2] < 0.) {
			   vect[2] = thick;
		   } else {
			   vect[2] = 0.;
		   }
			   
		   vect[0] = xsize * (rvec[0] - 0.5) + (vect[2] - thick/2.)*locdir[0]/locdir[2];
		   vect[1] = ysize * (rvec[1] - 0.5) + (vect[2] - thick/2.)*locdir[1]/locdir[2];
		   vect[3] = locdir[0]*pimom;
		   vect[4] = locdir[1]*pimom;
		   vect[5] = locdir[2]*pimom;
			   
/*  Propagate the track and make e-h pairs */

           deposit(vect, thick, new_drde, NEHSTORE, xeh, &neh);

/*  don't process overflows */

           if(neh < NEHSTORE) {

/*  Propagate the e's and h's */

              propag(thick, xsize, ysize, impwdth, implenn, implenp, temp, flux, rhe, rhh, ehole, neh, xeh, indeh);

/*  Count e's and h's on various pixels */

              detect(xsize, ysize, thick, impwdth, implenn, implenp, ehole, xeh, neh, pixel);

/*  Write out the results to a file */

              ofp = fopen(outfile, "a");
              fprintf(ofp,
              "%f %f %f %f %f %f %d \n", 
              vect[0], vect[1], vect[2], vect[3], vect[4], vect[5], neh);
              for (k=1; k < 2; ++k) {
                 for (j = 0; j < TYSIZE; ++j) {
                   fprintf(ofp,
                   "%2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f\n", 
                   pixel[k][0][j], pixel[k][1][j], pixel[k][2][j], pixel[k][3][j], pixel[k][4][j], 
			       pixel[k][5][j], pixel[k][6][j], pixel[k][7][j], pixel[k][8][j], pixel[k][9][j],
                   pixel[k][10][j], pixel[k][11][j], pixel[k][12][j], pixel[k][13][j], pixel[k][14][j],
                   pixel[k][15][j], pixel[k][16][j], pixel[k][17][j], pixel[k][18][j], pixel[k][19][j], 
			       pixel[k][20][j]);
                 }
              }
              fclose(ofp);          
		   } 
    
		   ievent += 1;
		  
/* Save current random number state */    
    
		   rluxut_(ivec);
		   isfp = fopen(seedfile, "w");
		   fprintf(isfp,
            "%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x%9x %d %d\n", 
            ivec[0], ivec[1], ivec[2], ivec[3], ivec[4], ivec[5], ivec[6], ivec[7],
            ivec[8], ivec[9], ivec[10], ivec[11], ivec[12], ivec[13], ivec[14], ivec[15],
            ivec[16], ivec[17], ivec[18], ivec[19], ivec[20], ivec[21], ivec[22], ivec[23],
            ivec[24], irun, ievent);
		   fclose(isfp);

/*  Determine current time */

		   now = time(NULL);
		   nows = localtime(&now);
		   sec = (*nows).tm_sec;
		   min = (*nows).tm_min;
		   hour = (*nows).tm_hour;
		   yday = (*nows).tm_yday;

		   if(ievent < 25) {printf("day %d at %02d:%02d:%02d, run %d, event %d, number of e-h pairs = %d\n", yday, hour, min, sec, irun, ievent, neh);}
		  
	   }
	  
	    irun += 1;
	    ievent = 0;
	}
   
/*  Determine current time */
   
   now = time(NULL);
   nows = localtime(&now);
   sec = (*nows).tm_sec;
   min = (*nows).tm_min;
   hour = (*nows).tm_hour;
   yday = (*nows).tm_yday;
   
   printf("End on day %d at %02d:%02d:%02d\n", yday, hour, min, sec);

    
} /* MAIN__ */ 


/* Subroutine */ int runinit(int *initseed, int *fileoff, int *runsize, float *lenxmin, float *lenxmax, float *deltaxlen, float *lenymin, float *lenymax, float *deltaylen)
{
    /* Initialize everything from external file */
    
    FILE *ifp;
	
	/* ******************************************************************** */
	/* * This routine initializes local run parameters from file run.init * */
	/* * Parameters: initseed - the random number seed                    * */
	/* *             fileoff - local offset for the first output file #   * */
	/* *             runsize - the # events in each run                   * */
	/* *             lenxmin - the cluster x-length of the first run      * */
	/* *                       (signed so that +len is along +x axis)     * */
	/* *             lenxmax - the cluster x-length of the last run       * */
	/* *           deltaxlen - the cluster x-length increment             * */
	/* *             lenymin - the cluster y-length of the first run      * */
	/* *                       (signed so that +len is along +y axis)     * */
	/* *             lenymax - the cluster y-length of the last run       * */
	/* *           deltaylen - the cluster y-length increment             * */
	/* ******************************************************************** */
	
    /* Function Body */
	
	/*  Initalize the parameters */
	
    ifp = fopen("run.init", "r");
    if (ifp==NULL) {
		printf("no run.init initialization file found/n");
		return 0;
    }
	
	/* next, the random number seed and cluster size parameters */    
    
    fscanf(ifp,"%d %d %d %f %f %f %f %f %f", initseed, fileoff, runsize, lenxmin, lenxmax, deltaxlen, lenymin, lenymax, deltaylen);
    
    printf("rand # seed = %d, file offset = %d, %d events/run, x-clust len from %f to %f in incr of %f, y-clust len from %f to %f in incr of %f\n", 
	       *initseed, *fileoff, *runsize, *lenxmin, *lenxmax, *deltaxlen, *lenymin, *lenymax, *deltaylen);
    
    fclose(ifp);
	
    return 0;
} /* runinit_ */

#include "pixelav3d.c"
