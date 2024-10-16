
/* *******************************************************************   
 * Original Author: Morris Swartz
 * Translator to C++: Oz Amram
 *
 *
 * This program reads and aanalyzes 21x13 pixel summary files.     *   
 * Create cluster length templates for x-/y-fitting procedure      *
 * Implement 2-pass appproach (3/6/06)                             *
 * Store templates at Nx x Ny                                       *
 * Assume that data are corrected for ROC response.                *
 * Use input thickness (zsize).                                    *
 * Add (0,0) point to each fit (8 Feb 2007).                       *
 * Fake z error parameters for short clusters if possible (8 Feb)  *
 * Reverse run processing direction to make faking work with cotb  *
 * reflection scheme (cotbeta > 0 only is stored temp) [1 Mar 07]  *
 * Add pixmax information to both output files.       [21 Oct 08]  *
 * Make adjustable array and template sizes.          [21 Oct 08]  *
 * Add digits to the error parameterization.          [22 Jul 09]  *
 * Modify small cluster error parameterization.       [22 Jul 09]  *
 * Add single pixel processing to estimate Lor widths [31 Oct 09]  *
 * Update for new response modeling input files       [28 Jan 10]  *
 *
 * Translated from fortran to C++                      May 2019
  * Increase size of buffers for Phase 2 and fix for Phase 1       
 *******************************************************************
 */
#include "template_utils.h"

void gen_xy_template2d(const int nevents = 32000, const int npt = 798, const int nprm = 5, 
        const int maxarg = 5){

    //size of the x and y arrays in template
    const int Nx = 21; // nbins in x
    const int Ny = 13; // nbins in y

    double xpar[nprm][4], xytemp[Nx][Ny][7][7];
    int nxytry[7][7];

    //declare larger arrays dynamically to avoid stack overflow
    float ***pixev = setup_3d_array(Nx, Ny, nevents); // charge per bin
    double *qsum = new double[nevents]; // sum charged
    int *i2d = new int[nevents]; // this is xhit position per bin
    int *j2d = new int[nevents]; // this is yhit position per bin

    //size of pixel (in pixelav coordinates)
    float xsize, ysize, zsize;

    // read pix_2t proc (containing # of template files to be processed)
    FILE *f_config = fopen("pix_2t.proc", "r");

    if(f_config ==0){
        printf("Can't find pix_2t.proc! Exiting \n");
        exit(1);
    }

    int file_start, num_files, non_linear;
    float q100, q101, q100_frac, noise, common_frac, gain_frac, readout_noise;
    fscanf(f_config,"%d %d %f %f %f %f %f %f %f %d", &file_start, &num_files, &noise, &q100, &q101, &q100_frac, &common_frac, 
            &gain_frac, &readout_noise, &non_linear);
    fclose(f_config);
    printf("processing %d files starting from %d, noise = %f, threshold0 = %f, threshold1 = %f," 
            "rms threshold frac = %f, common_frac = %f, gain fraction = %f, readout noise = %f, nonlinear_resp = %d \n", 
            num_files, file_start, noise, q100, q101, q100_frac, common_frac, gain_frac, readout_noise, non_linear);


    //define RMS of noise
    double rten = 10.;
    double thr = q100;
    double thr10 = 0.1*thr;

    char fname[100];
    char header[120];
    char buffer[120];

    //loop over files backwards
    for(int iFile = file_start + num_files - 1; iFile >= file_start; iFile--){
        sprintf(fname, "template_events_d%05i.out", iFile);
        FILE *f_evts = fopen(fname, "r");

        if(f_evts ==0) {
            printf("Cant find file %s! Exiting \n", fname);
            exit(1);
        }

        fgets(header, 120, f_evts);

        fgets(buffer, 120, f_evts); 

        sscanf(buffer, "%f %f %f", &xsize, &ysize, &zsize);
        printf("xsize = %f,  ysize = %f, zsize = %f \n", xsize, ysize, zsize);

        float z_center = zsize/2.0;

        float cosx, cosy, cosz, x1,y1,z1;

        double qavg=0.;

        int nelec;

        zero_3d_array(pixev, Nx, Ny, nevents);

        //loop over events 
        int n;
        for(n=0; n < nevents; n++){
            if(fgets(buffer, 120, f_evts) ==0)
                //end of file
                break;

            sscanf(buffer, "%f %f %f %f %f %f %i",
                    &x1, &y1, &z1, &cosx, &cosy, &cosz, &nelec);

            float bixin[Nx][Ny];
            memset(bixin, 0., sizeof(bixin));

            //note the confusing index switch. input arrays are 13x21, we swap
            //to 21x13
            for(int j=0; j<Ny; j++){ // 13
                for (int i=0; i<Nx; i++){ // 21
                    fscanf(f_evts, " %f ", &bixin[i][j]); // this is charge in 1000s of e-
                }
            }

            //  Pixelav gives hit position at face of pixel, translate to
            //  3d center of the pixel
            //  propagate from edge to center of pixel
            float xhit = x1 + (z_center - z1) * cosx/cosz;
            float yhit = y1 + (z_center - z1) * cosy/cosz;

            //hit position bins (7 bins within a pixel)
            //last bin is the same as first, translated 1 pixel over
            // i.e. bin0 px0 = bin7 px-1 
            // and  bin7 px0 = bin0 px1
            //divide into 2d bins (6x6) for max charge

            // restrict from 0 to 6;    
            int i2dn  = int(xhit/xsize*6 + 3.5);
            i2dn = std::min(6, std::max(i2dn, 0));
            i2d[n] = i2dn;

            //restrict from 0 to 6; 
            int j2dn  = int(yhit/ysize*6 + 3.5);
            j2dn = std::min(6, std::max(j2dn, 0));
            j2d[n] = j2dn;

            qsum[n] = 0.;

            for(int i=0; i<Nx; i++){
                for(int j=0; j<Ny; j++){
                    double qin = bixin[i][j] * rten;
                    if(qin < 0.) qin = 0.;
                    pixev[i][j][n] = qin;
                    qsum[n] += qin;
                }
            }
            qavg += qsum[n];

        } // finish loop over events

        int nmc = n; // total number of events

        printf("finish reading run %i read %i events, qsum = %.1f \n", iFile, nmc, qavg);

        fclose(f_evts);
        qavg /= double(nmc); // get average charge
        printf("finish reading run %i read %i events, qavg = %.1f \n", iFile, nmc, qavg);

        //zero things for this run
        memset(nxytry, 0, sizeof(nxytry));
        memset(xytemp, 0., sizeof(xytemp));

        for(n=0; n<nmc; n++){

            //count number of pixels with less than average charge
            //We make templates with only events where charge < average charge
            //to avoid long tails (charge is truncated in reco anyway)
            if(qsum[n] < qavg){

                //fill 2d 7x7 2d template for average charge determination
                //hit position indexes
                int k = i2d[n]; 
                int l = j2d[n]; 
                nxytry[k][l] +=1;

                for(int i=0; i < Nx; i++){
                    for(int j=0; j < Ny; j++){
                        xytemp[i][j][k][l] += pixev[i][j][n];
                    }
                }

                for(int i=0; i < Nx; i++){
                    for(int j=0; j < Ny; j++){
                        // 0th bin of new pixel same as last bin of previous pixel   
                        if(k==0 && i>0){ // bin0px(i) = bin6px(i-1) 
                            xytemp[i][j][6][l] += pixev[i-1][j][n]; // adding to bin6 of that same pix 
                        }
                        else if(k==6 && i>0){ // bin6px(i) = bin0px(i+1) 
                            xytemp[i-1][j][0][l] += pixev[i][j][n];  // adding to bin0 of pix-1 
                        }

                        if(l==0 && j>0){
                            xytemp[i][j][k][6] += pixev[i][j-1][n]; // adding to bin6 of that same pix 
                        }
                        else if(l==6 && j>0){
                            xytemp[i][j-1][k][0] += pixev[i][j][n]; // adding to bin0 of pix-1
                        }

                        if(i>0 && j>0) {
                            if(k==0 && l==0){
                                xytemp[i][j][6][6] += pixev[i-1][j-1][n];
                            }
                            else if(k==6 && l==6){
                                xytemp[i-1][j-1][0][0] += pixev[i][j][n];
                            }
                            else if(k==0 && l==6){
                                xytemp[i][j-1][6][0] += pixev[i-1][j][n];
                            }
                            else if(k==6 && l==0){
                                xytemp[i-1][j][0][6] += pixev[i][j-1][n];
                            }
                        }

                    } // end Ny
                } // end Nx

                if(k < 0 || k > 6){
                    printf("Problem with event %i, index (k) is =%i \n", n, k);
                }
                else if(k==0) nxytry[6][l] += 1;
                else if(k==6) nxytry[0][l] += 1;

                if(l < 0 || l > 6){
                    printf("Problem with event %i, index (l) is =%i \n", n, l);
                }
                else if(l==0) nxytry[k][6] += 1;
                else if(l==6) nxytry[k][0] += 1;

                if(k==0 && l==0){
                    nxytry[6][6] += 1;
                }
                else if(k==6 && l==6){
                    nxytry[0][0] += 1;
                }
                else if(k==0 && l==6){
                    nxytry[6][0] += 1;
                }
                else if(k==6 && l==0){
                    nxytry[0][6] += 1;
                }
            } 
        }//end loop over events

        //renorm the distributions, find the maximum average pixel signal in 7x7 grid
        float pixmax = 0.;
        int imin = Nx;
        int imax = 0;
        int jmin = Ny;
        int jmax = 0;
        for(int k =0; k<7; k++){
            for(int l =0; l<7; l++){
                for(int i =0; i<Nx; i++){
                    for(int j =0; j<Ny; j++){
                        float sigxy = xytemp[i][j][k][l] / float(nxytry[k][l]);
                        xytemp[i][j][k][l] = sigxy;
                        if(sigxy > pixmax) pixmax = sigxy;
                        if(sigxy > thr10) {
                            if(i < imin) imin=i; 
                            if(i > imax) imax=i;
                            if(j < jmin) jmin=j;
                            if(j > jmax) jmax=j;
                        }
                    }
                }
            }
        }
        //increase by 1 to match previous fortran indexing, this is accounted
        //for in the reco code
        imin = imin+1;
        imax = imax+1;
        jmin = jmin+1;
        jmax = jmax+1;

        // set fit parameters to 0 because they are unused 
        for(int p=0; p<nprm; p++){
            for(int i=0; i<4; i++){
                xpar[p][i] = 0;
            }
        }

        //open summary files
        char file_out[100];

        //x template file
        sprintf(file_out, "zptemp_%05i.txt", iFile);
        FILE *zptemp_file = fopen(file_out, "w+");
        fprintf(zptemp_file, "%9.6f %9.6f %9.6f \n", cosx, cosy, cosz);
        fprintf(zptemp_file, "%8.1f %8.1f %2i %2i %2i %2i\n", qavg,pixmax,imin,imax,jmin,jmax); 
        for(int k=2; k <4; k++){
            for(int p=0; p<nprm; p++){
                fprintf(zptemp_file, "%15.8E ", xpar[p][k]); // these should be 0 now
            }
            fprintf(zptemp_file, "\n");
        }

        for(int l=0; l <7; l++){
            for(int k=0; k <7; k++){
                float ycenter = ((l+1)*0.166667 - 0.666667) *ysize;
                float xcenter = ((k+1)*0.166667 - 0.666667) *xsize;
                fprintf(zptemp_file, "biny %2i,  ycenter = %8.2f um, binx %2i, xcenter = %8.2f um \n", l+1, ycenter, k+1, xcenter);
                for(int j=0; j<Ny; j++){
                    for(int i=0; i<Nx; i++){
                        fprintf(zptemp_file, "%8.1f ", xytemp[i][j][k][l] );
                    }
                    fprintf(zptemp_file, "\n"); 
                }
            }
        }

        fclose(zptemp_file);

    } //end loop over files

    delete_3d_array(pixev, Nx, Ny, nevents);
    delete[] qsum; 
    delete[] i2d; 
    delete[] j2d;
}


int main(){
    gen_xy_template2d();
    return 0;
}

