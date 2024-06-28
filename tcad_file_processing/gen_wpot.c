/* Program to interpolate regular efield grid from irregular mesh in TCAD files */
/* Usage: gen_wpot_rev_plt tcad_dir [for example gen_efield_rev_plt dot1_150x100_prod_dj44] */
/* inputs are x,y,z grid node numbers desired and option to zero efield starting at some z (use 0 for this) */
/* Does 3-D interpolation using nearest 4 mesh vertices that are non-planar (can fail requiring NMMAX and NMMIN to be increased) */
/*also transforms from TCAD coordinate to pixelav coordinates */
/* Output: file weighting.out contains the map needed to make wgt_pot.init files for pixelav */
/*         file wplot.out contains a profile along the symmetry axis (center of cell) with mesh point spacing for diagnostic purpose */
/* adjust the zeroed section to simulate he low field region of a p+ in n strip detector */

#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#define SWAP(a,b)  temp=(a); (a)=(b); (b)=temp;
#define M 7
#define NSTACK 50 
#define NMMAX 3000 /*number of nearest neighbors to search for 4 non-planar nearest neighbors in busy part of lattice (near n+ side) */
#define NMMIN 600 /*number of nearest neighbors to search for 4 non-planar nearest neighbors in unstructured part of sensor (near p+ side) */
#define NMVTX 150000 /*max number of vertices */
int main(int argc, char *argv[])
{
    static int i, j, k, l, m, n, ia, nvtx, nevals, nx, ix[NMMAX], ixcop[NMMAX], ixtmp;
    static float xvtx[3][NMVTX], wvtx[NMVTX], xaxis, yaxis, tol, halfnum;
    static char gridfile[100], desfile[100], inp1[100], inp2[100], inp3[100];
    FILE *gridfp, *desfp, *effp, *pltfp;
    static double size[3] = {20.,20.,100.}, ds[3], x[3], wgt, dxtmp, delta, dotp1, dotp2;
    double xp[2], xarray[9][2], wgtarray[9], xc, yc, xt, yt, xr, yr;
    static double a[3][3], ainv[3][3], df[3], alpha[3], det, atest[3][3];
    static double x01[3], x02[3], x03[3], cross[3], a01, a02, a03, ac;
    static double xmin[3], xmax[3];
    float dx[NMMAX], ixs[NMMAX], dxcop[NMMAX];
    float *dxp, *ixsp;
    static int ns[3] = {21,21,92}, ninp[3];
    static unsigned long nsort;
    static double temp=283.;
    static double zold;

/* Check number of arguments */

    if(argc != 3) {
      printf("wrong number of arguments = %d\n",argc);
      return 1;
    }

/* initialize pointers for the quicksort */
    
    dxp = dx - 1;
    ixsp = ixs - 1;
    
   /* Construct file names */
   strcpy(gridfile, argv[1]);
   strcat(gridfile,"/");
   strcat(gridfile, argv[1]);
   strcat(gridfile,"_msh.grd");
   strcpy(desfile, argv[1]);
   strcat(desfile,"/");
   strcat(desfile, argv[1]);
   strcat(desfile,"_");
   strcat(desfile, argv[2]);
   strcat(desfile,"_des.dat");
   printf("Grid file = %s, dessis plot file = %s\n", gridfile, desfile);
   
   /* Make sure files are available */
   
   gridfp = fopen(gridfile, "r");
   if (gridfp==NULL) {
      printf("can't find %s\n",gridfile);
      return 1;
   }
    
/* find and read in the grid coordinates */

       for(m=0;m<3;++m){
          xmin[m]=0.;
          xmax[m]=0.;
       }
       while (fscanf(gridfp,"%s", inp1) !=EOF) {
         if(strcmp(inp1,"Vertices") == 0) {
           fscanf(gridfp,"%s %d %s %s", inp1, &nvtx,inp2,inp3);
           printf("number of vertices = %d\n",nvtx);
           if(nvtx > NMVTX) {printf("too many vertices \n"); return 0;}
           for(i=0; i<nvtx; ++i) {
              fscanf(gridfp,"%f %f %f", &xvtx[2][i], &xvtx[0][i], &xvtx[1][i]);
              for(m=0;m<3;++m){
                 if(xvtx[m][i] < xmin[m]) {xmin[m] = xvtx[m][i];}
                 if(xvtx[m][i] > xmax[m]) {xmax[m] = xvtx[m][i];}
              }
           } 
           goto closeg;
          }
        }
        
/* close grid file */
        
closeg: fclose(gridfp);
    
    desfp = fopen(desfile, "r");
    if (desfp==NULL) {
       printf("can't find %s\n",desfile);
       return 1;
    }
    
/* find and read in the potential at each vertex */

       while (fscanf(desfp,"%s", inp1) !=EOF) {
         if(strcmp(inp1,"(\"ElectrostaticPotential\")") == 0) {
           while (fscanf(desfp,"%s", inp1) !=EOF) {
             if(strcmp(inp1,"Values") == 0) {
               while (fscanf(desfp,"%s", inp1) !=EOF) {
                 if(strcmp(inp1,"{") == 0) {
                   for(i=0; i<nvtx; ++i) {
                      fscanf(desfp,"%f", &wvtx[i]);
                   } 
                   goto closed;
                }
              }
             }
           }
          }
        }
        
/* close dessis plot file */
        
closed: fclose(desfp);
   printf("first node = %f, last node = %f \n", wvtx[0], wvtx[nvtx]);
        
/* open the efield output file */
        
       effp = fopen("weighting.out", "w");
        
       pltfp = fopen("wplot.out", "w");
       
/* calculate the size of the detector and determine the number of output grid points */
            
               
       for(l=0; l<3; ++l) {
            size[l] = xmax[l] - xmin[l];
       }       
       printf("detector dimensions = %f %f %f um \n", size[0], size[1], size[2]);
       printf("enter the number of 1/2 pixels in each direction [5]\n");
       scanf("%f", &halfnum);
       
/* calculate the size of the detector and determine the number of output grid points */
               
       for(l=0; l<2; ++l) {
          size[l] /= halfnum;
       }
       
       printf("enter the number of output grid points in each dim: nx[21], ny[21], nz[92]\n");
       scanf("%d %d %d", &ninp[0], &ninp[1], &ninp[2]);
       for(m=0;m<3;++m){
         if(ninp[m]>0) {ns[m] = ninp[m];}
       }
       printf("(nx, ny, nz) = (%d, %d, %d)\n", ns[0], ns[1], ns[2]);
       
       printf("enter the x y coordinates and tolerance for plot axis\n");
       scanf("%f %f %f", &xaxis, &yaxis, &tol);

/* Decide if we need to zero some of the p-side to avoid double junction problems (in inverted Si) */

       l=0;
       for(i=0; i<nvtx; ++i) {
          if(fabs(xvtx[0][i]-xaxis) < tol && fabs(xvtx[1][i]-yaxis) < tol) {
            ixs[l]= (float) i;
            dx[l]=xvtx[2][i];
            ++l;
            if(l>=NMMAX) {break;}
           }
        }
        nsort = l-1;
        if(nsort <= 0) {
           printf("nsort = %ld, no central elements found\n", nsort);
           return 1;
        }
   
   printf("nsort = %ld \n", nsort);

/* sort indices into increasing z */

        sort2(nsort, dxp, ixsp);
      
/* print increasing z points and the electric field */

        for(m=0; m<nsort; ++m) {
            ixtmp = (int) ixs[m];
            printf("%d ind = %d, z = %e, wvtx = %e\n", m, ixtmp, xvtx[2][ixtmp], wvtx[ixtmp]);
            fprintf(pltfp,"%e %e\n", xmax[2]-xvtx[2][ixtmp], wvtx[ixtmp]);
           zold = xvtx[2][ixtmp];
        }
       fclose(pltfp);


/* Now create a regular grid of interpolated values for use in pixelav */

         for(l=0; l<3; ++l) {
            ds[l] = size[l]/((float) (ns[l]-1));
         }
         
         for(k=0; k<ns[2]; ++k) {
            x[2] = xmin[2]+k*ds[2];
            
            if(x[2]>(size[2]-30.) || x[2]<30.) {
               nsort = NMMAX;
            } else {
               nsort = NMMIN;
            }
            for(j=0; j<ns[1]; ++j) {
               xp[1] = j*ds[1];
               for(i=0; i<ns[0]; ++i) {
                 xp[0] = i*ds[0];
                 
/* Generate 9 sets of points to interpolate */

                 xc = xp[0]; xt = xp[0]+2.*size[0]; xr = 2.*size[0] - xp[0];
                 yc = xp[1]; yt = xp[1]+2.*size[1]; yr = 2.*size[1] - xp[1];
                 xarray[0][0] = xt; xarray[0][1] = yr; 
                 xarray[1][0] = xc; xarray[1][1] = yr; 
                 xarray[2][0] = xr; xarray[2][1] = yr; 
                 xarray[3][0] = xt; xarray[3][1] = yc; 
                 xarray[4][0] = xc; xarray[4][1] = yc;
                 xarray[5][0] = xr; xarray[5][1] = yc; 
                 xarray[6][0] = xt; xarray[6][1] = yt; 
                 xarray[7][0] = xc; xarray[7][1] = yt; 
                 xarray[8][0] = xr; xarray[8][1] = yt;
                  
                 for(ia=0; ia < 9; ++ia) {
                    x[0] = xarray[ia][0];
                    x[1] = xarray[ia][1];
                 
/* Search for the 4 nearest non-planar vertices, first load NMEM differences */

                    for(l=0; l<nsort; ++l) {
                       delta = sqrt((double)((x[1]-xvtx[1][l])*(x[1]-xvtx[1][l]) + (x[2]-xvtx[2][l])*(x[2]-xvtx[2][l])
                                 +(x[0]-xvtx[0][l])*(x[0]-xvtx[0][l])));
                       ixs[l]= (float) l;
                       dx[l]=(float) delta;
                    }
                 
/* Next, sort them */

                   sort2(nsort, dxp, ixsp);
                   for(m=0; m<nsort; ++m) {
                      ix[m] = (int) ixs[m];
                   }
                 
/* Next, continue through the list */

                     for(l=nsort; l<nvtx; ++l) {
                        delta = sqrt((double)((x[1]-xvtx[1][l])*(x[1]-xvtx[1][l]) + (x[2]-xvtx[2][l])*(x[2]-xvtx[2][l])
                                 +(x[0]-xvtx[0][l])*(x[0]-xvtx[0][l])));
                        if(delta>dx[nsort-1]) {continue;}
                        for(m=0;m<nsort;++m) {             
                           if(delta<dx[m]) {
                              if(m<(nsort-1)) {
                                 for(n=(nsort-1); n>m; --n) {
                                    dx[n] = dx[n-1];
                                    ix[n] = ix[n-1];
                                 }
                              }
                              dx[m] = (float) delta;
                              ix[m] = l;
                              break;  
                           }
                        }
                     }
                 
/* save the sorted list for later diagnotic use */
                 
                    for(m=0; m<nsort; ++m) {
                       ixcop[m] = ix[m];
                       dxcop[m] = dx[m];
                    }

                 
/* nx is the size of the closest neighbor stack */                 
                 
                    nx=nsort;
                 
/* search the list for non-planar vectors. first create a vector from the closest pair */                 
                 
                   for(m=0;m<3;++m){ x01[m] = (double)(xvtx[m][ix[1]] - xvtx[m][ix[0]]); }
                   a01 = sqrt(x01[0]*x01[0]+x01[1]*x01[1]+x01[2]*x01[2]);  
                  
/* look for a vector that is not co-linear */  
               
secondvec:         for(m=0;m<3;++m){ x02[m] = (double)(xvtx[m][ix[2]] - xvtx[m][ix[0]]); }
                   a02 = sqrt(x02[0]*x02[0]+x02[1]*x02[1]+x02[2]*x02[2]);
                   dotp1=0.;
                   for(m=0;m<3;++m){ dotp1 += x02[m]*x01[m]; }
                   dotp1=dotp1/(a01*a02);
                   if(fabs(dotp1) > 0.98) {
                 
/* drop the stack by one */

                     if (nx <= 4) {
                         printf(" no 2nd vec neighbors, nx = %d\n", nx);
                         for(m=0;m<nsort;++m) {
                            printf("%d %d dx = %e, xvtx = %e %e %e \n", m, ixcop[m], dxcop[m], 
                               xvtx[0][ixcop[m]], xvtx[1][ixcop[m]], xvtx[2][ixcop[m]]);
                         }
                         return 1;
                     } else {
                         for(m=2; m<(nx-1); ++m) {
                            dx[m] = dx[m+1];
                            ix[m] = ix[m+1];
                         }
                         nx=nx-1;
                         goto secondvec;
                      }
                   }
                  
/* look for a vector that is not co-planar */ 

                   cross[0] = x01[1]*x02[2]-x01[2]*x02[1];
                   cross[1] = x01[2]*x02[0]-x01[0]*x02[2];                
                   cross[2] = x01[0]*x02[1]-x01[1]*x02[0];                
                   ac = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]); 
thirdvec:          for(m=0;m<3;++m){ x03[m] = (double)(xvtx[m][ix[3]] - xvtx[m][ix[0]]); }
                   a03 = sqrt(x03[0]*x03[0]+x03[1]*x03[1]+x03[2]*x03[2]);
                   dotp2=0.;
                   for(m=0;m<3;++m){ dotp2 += x03[m]*cross[m]; }
                   dotp2=dotp2/(ac*a03);
                   if(fabs(dotp2) < 0.01 ) {
                 
/* drop the stack by one */

                      if (nx <= 4) {
                         printf(" no 3rd vec neighbors, nx = %d\n", nx);
						 printf(" grid coords = %e %e %e \n", x[0], x[1], x[2]);
                         printf("a01 = %e, x01 = %e %e %e \n",a01, x01[0], x01[1], x01[2]);
                         printf("a02 = %e, x01 = %e %e %e \n",a02, x02[0], x02[1], x02[2]);
                         printf("a03 = %e, x01 = %e %e %e \n",a03, x03[0], x03[1], x03[2]);
                         printf("ac = %e, cross = %e %e %e \n",ac, cross[0], cross[1], cross[2]);
                         printf("dotp1 = %e, dotp2 = %e\n", dotp1, dotp2);
                         for(m=0;m<nsort;++m) {
                            printf("%d %d dx = %e, xvtx = %e %e %e \n", m, ixcop[m], dxcop[m], 
                               xvtx[0][ixcop[m]], xvtx[1][ixcop[m]], xvtx[2][ixcop[m]]);
                         }
                         return 1;
                      } else {
                         for(m=3; m<(nx-1); ++m) {
                            dx[m] = dx[m+1];
                            ix[m] = ix[m+1];
                      }
                      nx=nx-1;
                      goto thirdvec;
                    }
                  }  
                 
/* We now have a list of the 4 closest vertices, let's interpolate */

                  for(m=0; m<3; ++m) {
                     a[0][m] = x01[m];
                     a[1][m] = x02[m];
                     a[2][m] = x03[m];
                  }
                 
/* Calculate the determinant and check to see that it's non-zero (the points are non-planar) */

                  det = a[0][0]*a[1][1]*a[2][2]+a[0][1]*a[1][2]*a[2][0]+a[0][2]*a[1][0]*a[2][1]
                     - a[2][0]*a[1][1]*a[0][2]-a[2][1]*a[1][2]*a[0][0]-a[2][2]*a[1][0]*a[0][1];
                     
                  if(fabs(det) < 1.e-6) {
                 
/* If all neighbors are planar, drop the stack by one */

                     printf(" determinant is very small = %f\n", det);
                     return 1;
                  }
                 
/* Calculate the inverse of the matrix a */

                  ainv[0][0] = (a[1][1]*a[2][2]-a[2][1]*a[1][2])/det;
                  ainv[1][0] =-(a[1][0]*a[2][2]-a[2][0]*a[1][2])/det;
                  ainv[2][0] = (a[1][0]*a[2][1]-a[2][0]*a[1][1])/det;
                  ainv[0][1] =-(a[0][1]*a[2][2]-a[2][1]*a[0][2])/det;
                  ainv[1][1] = (a[0][0]*a[2][2]-a[2][0]*a[0][2])/det;
                  ainv[2][1] =-(a[0][0]*a[2][1]-a[2][0]*a[0][1])/det;
                  ainv[0][2] = (a[0][1]*a[1][2]-a[1][1]*a[0][2])/det;
                  ainv[1][2] =-(a[0][0]*a[1][2]-a[1][0]*a[0][2])/det;
                  ainv[2][2] = (a[0][0]*a[1][1]-a[1][0]*a[0][1])/det;
              
/* Interpolate each component of the electric field */

                  df[0] = (double) (wvtx[ix[1]] - wvtx[ix[0]]);
                  df[1] = (double) (wvtx[ix[2]] - wvtx[ix[0]]);
                  df[2] = (double) (wvtx[ix[3]] - wvtx[ix[0]]);
                  for(m=0;m<3;++m){
                     alpha[m] = 0.;
                     for(n=0;n<3;++n){ alpha[m] += ainv[m][n]*df[n];}
                  }
                  wgt=wvtx[ix[0]];
                  for(n=0;n<3;++n){ wgt += alpha[n]*(x[n]- (double)xvtx[n][ix[0]]);}
                  wgtarray[ia] = wgt;
                } 
                fprintf(effp,"%4d%4d%4d %e %e %e %e %e %e %e %e %e\n", i+1, j+1, k+1, 
                wgtarray[0],wgtarray[1],wgtarray[2],wgtarray[3],wgtarray[4],wgtarray[5],
                wgtarray[6],wgtarray[7],wgtarray[8]);
              }
            }
        }
        return 0;
}
