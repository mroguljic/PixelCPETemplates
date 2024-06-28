/* Program to interpolate regular efield grid from irregular mesh in TCAD files */
/* Usage: gen_efield_rev_plt tcad_dir voltage [for example gen_efield_rev_plt dot1_150x100_prod_dj44 300] */
/* inputs are x,y,z grid node numbers desired and option to zero efield starting at some z (use 0 for this) */
/* Does 3-D interpolation using nearest 4 mesh vertices that are non-planar (can fail requiring NMMAX and NMMIN to be increased) */
/*also transforms from TCAD coordinate to pixelav coordinates */
/* Output: file efield.out contains the map needed to make pixel.init files for pixelav */
/*         file eplot.out contains an E_z profile along the symmetry axis (center of cell) with mesh point spacing for diagnostic purpose */
/* adjust the zeroed section to simulate he low field region of a p+ in n strip detector */
/* Version for 2-fold symmetry instead of 4-fold symmetry */

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
int main(int argc, char *argv[])
{
    static int i, j, k, l, m, n, nvtx, nevals, nx, ix[NMMAX], ixcop[NMMAX], ixtmp;
    static float xvtx[3][30000], evtx[3][30000];
    static char gridfile[100], desfile[100], inp1[100], inp2[100], inp3[100];
    FILE *gridfp, *desfp, *effp, *pltfp;
    static double size[3] = {20.,20.,100.}, ds[3], x[3], e[3], dxtmp, delta, dotp1, dotp2;
    static double a[3][3], ainv[3][3], df[3], alpha[3], det, atest[3][3];
    static double x01[3], x02[3], x03[3], cross[3], a01, a02, a03, ac;
    static double xmin[3], xmax[3];
    float dx[NMMAX], ixs[NMMAX], dxcop[NMMAX], zmin;
    float *dxp, *ixsp;
    static int ns[3] = {21,21,92}, ninp[3];
    static unsigned long nsort;
    static double temp=283.;
    static double mu, vm, ec, beta, ibeta, xsum, esum, zsum, zold;

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
    
/* find and read in the electric field at each vertex */

       while (fscanf(desfp,"%s", inp1) !=EOF) {
         if(strcmp(inp1,"(\"ElectricField-Vector\")") == 0) {
           while (fscanf(desfp,"%s", inp1) !=EOF) {
             if(strcmp(inp1,"Values") == 0) {
               while (fscanf(desfp,"%s", inp1) !=EOF) {
                 if(strcmp(inp1,"{") == 0) {
                   for(i=0; i<nvtx; ++i) {
                      fscanf(desfp,"%f %f %f", &evtx[2][i], &evtx[0][i], &evtx[1][i]);
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
        
/* open the efield output file */
        
       effp = fopen("efield.out", "w");
        
       pltfp = fopen("eplot.out", "w");
       
/* calculate the size of the detector and determine the number of output grid points */
               
       for(l=0; l<3; ++l) {
            size[l] = xmax[l] - xmin[l];
       }
       
       printf("detector dimensions = %f %f %f um \n", size[0], size[1], size[2]);
       printf("enter the number of output grid points in each dim: nx[21], ny[21], nz[92]\n");
       scanf("%d %d %d", &ninp[0], &ninp[1], &ninp[2]);
       for(m=0;m<3;++m){
         if(ninp[m]>0) {ns[m] = ninp[m];}
       }
       printf("(nx, ny, nz) = (%d, %d, %d)\n", ns[0], ns[1], ns[2]);

/* Decide if we need to zero some of the p-side to avoid double junction problems (in inverted Si) */

       l=0;
       for(i=0; i<nvtx; ++i) {
          if(xvtx[0][i]==0. && xvtx[1][i]==(xmin[1]+size[1]/2.)) {
            ixs[l]= (float) i;
            dx[l]=xvtx[2][i];
            ++l;
            if(l>=NMMAX) {break;}
           }
        }
        nsort = l-1;
        if(nsort <=0) {
           printf("nsort = %ld, no central elements found\n", nsort);
           return 1;
        }

/* sort indices into increasing z */

        sort2(nsort, dxp, ixsp);
      
   /* calculate mobility constants */
   
   vm=1.53e9 * pow(temp, -0.87);
   ec=1.01 * pow(temp, 1.55);
   beta = 2.57e-2 * pow(temp, 0.66);
   ibeta=1./beta;
   xsum = 0.;
   esum = 0.;
   zsum = 0.;

/* print increasing z points and the electric field */

        for(m=0; m<nsort; ++m) {
            ixtmp = (int) ixs[m];
            printf("%d ind = %d, z = %e, evtx = %e %e %e\n", m, ixtmp, xvtx[2][ixtmp], evtx[0][ixtmp], 
                    evtx[1][ixtmp], evtx[2][ixtmp]);
            fprintf(pltfp,"%e %e\n", xmax[2]-xvtx[2][ixtmp], fabs(evtx[2][ixtmp]));
           if(m > 0) {
              mu=vm/ec/pow((1.+pow(fabs(evtx[2][ixtmp])/ec,beta)), ibeta);
              xsum += 1.086*mu*3.8e-4*(xvtx[2][ixtmp] - zold);
              esum += evtx[2][ixtmp]*(xvtx[2][ixtmp] - zold);
              zsum += (xvtx[2][ixtmp] - zold);
           }
           zold = xvtx[2][ixtmp];
        }
       fclose(pltfp);
       printf("total Lorentz Drift = %lf, average efield = %lf \n", xsum, esum/zsum);
       printf("enter zmin (E = 0 for z>zmin) \n");
       scanf("%f", &zmin);
       printf("zmin = %f um\n",zmin);        

/* Now create a regular grid of interpolated values for use in pixelav */

         for(l=0; l<3; ++l) {
            ds[l] = size[l]/((float) (ns[l]-1));
         }
         
         for(k=0; k<ns[2]; ++k) {
            x[2] = xmin[2]+k*ds[2];
            
/* offset the end points to get non-zero field */

            if(k==0) {x[2] = x[2] + 1.2;}
            if(k==(ns[2]-1)) {x[2] = x[2] - 1.2;}
            if(x[2]>(size[2]-30.) || x[2]<30.) {
               nsort = NMMAX;
            } else {
               nsort = NMMIN;
            }
            for(j=0; j<ns[1]; ++j) {
               x[1] = xmin[1]+j*ds[1];
               for(i=0; i<ns[0]; ++i) {
                 x[0] = xmin[0]+i*ds[0];
                 if(x[2] < zmin) {
                 
/* Search for the 4 nearest non-planar vertices, first load NMEM differences */

                 for(l=0; l<nsort; ++l) {
                    delta = sqrt((double)((x[1]-xvtx[1][l])*(x[1]-xvtx[1][l]) + (x[2]-xvtx[2][l])*(x[2]-xvtx[2][l])
                                 +(x[0]-xvtx[0][l])*(x[0]-xvtx[0][l])));
                    ixs[l]= (float) l;
                    dx[l]=(float) delta;
                 }
                 
/* Next, sort them */
/*  bsort:           dxtmp=0.;
                 for(m=0; m<49; ++m){
                    if(dx[m] > dx[m+1]) {
                    dxtmp=dx[m];
                    ixtmp=ix[m];
                    dx[m]=dx[m+1];
                    ix[m]=ix[m+1];       
                    dx[m+1]=dxtmp;
                    ix[m+1]=ixtmp;
                    goto bsort;
                  }
                }
 */
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
               
secondvec:     for(m=0;m<3;++m){ x02[m] = (double)(xvtx[m][ix[2]] - xvtx[m][ix[0]]); }
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
thirdvec:       for(m=0;m<3;++m){ x03[m] = (double)(xvtx[m][ix[3]] - xvtx[m][ix[0]]); }
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

                 for(l=0; l<3; ++l) {
                    df[0] = (double) (evtx[l][ix[1]] - evtx[l][ix[0]]);
                    df[1] = (double) (evtx[l][ix[2]] - evtx[l][ix[0]]);
                    df[2] = (double) (evtx[l][ix[3]] - evtx[l][ix[0]]);
                    for(m=0;m<3;++m){
                       alpha[m] = 0.;
                       for(n=0;n<3;++n){ alpha[m] += ainv[m][n]*df[n];}
                    }
                    e[l]=evtx[l][ix[0]];
                    for(n=0;n<3;++n){ e[l] += alpha[n]*(x[n]- (double)xvtx[n][ix[0]]);}
                 }
                } else {
                  e[0] = 0.; e[1] = 0.; e[2] = 0.;
                }
                fprintf(effp,"%4d%4d%4d   %e %e %e \n", i+1, j+1, k+1, e[0],e[1],e[2]);
              }
            }
        }
        return 0;
}
