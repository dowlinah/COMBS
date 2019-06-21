/* --- appbt.c ---
 *
 * NAS CFD code: APPBT
 *
 * Original FORTRAN version:
 *
 * Author: Sisira Weeratunga
 *         NASA Ames Research Center
 *         (10/25/90)
 *
 * This is the C version of the FORTRAN code developed at NASA Ames Research
 * Center.  We converted the serial code to C as a part of our class project
 * for CS 838-3/ChE 562, which was offered in Spring 1993 by Mark Hill, 
 * Sangtae Kim, and Mary Vernon at the University of Wisconsin-Madison.  
 *
 *					July 20, 1993.
 *
 *					Doug Burger
 *					dburger@mysost.cs.wisc.edu
 *					Computer Sciences Department
 *
 *					Sanjay Mehta
 *					sanjaym@luther.che.wisc.edu
 *					Chemical Engineering Department
 *					
 *					University of Wisconsin-Madison
 *
 */
  
#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <math.h>
#include "appbt.h"

Global_Struct *GMEM;

int nx, ny, nz;
int ipr, inorm, itmax;
int subblocks[3];
FILE *in_file, *out_file;
clock_t looptotal, alltotal; /* Timing variables */
double rsdnm[5];
double dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, dy5, dz1, dz2, dz3, dz4, dz5;

int setbv_()
{
  int i, j, k;
  extern int exact_();
  int x0, x1, xn, y0, y1, yn, z0, z1, zn;

/* ***set the boundary values of dependent variables */
  
/* ***set the dependent variable values along the top and bottom faces */
  
  XO(x1, 1)
  YO(y1, 1)
  ZO(z1, 1)
  XO(xn, nx)
  YO(yn, ny)
  ZO(zn, nz)
  
  for (j = 1; j <= ny; ++j) {
    YO(y0, j)   
    for (i = 1; i <= nx; ++i) {
      XO(x0, i)      
      exact_(i, j, 1, &UA[5 * (x0 + y0 + z1)]);
      exact_(i, j, nz, &UA[5 * (x0 + y0 + zn)]);
    }
  }
  
/* ***set the dependent variable values along north and south faces */
  
  for (k = 1; k <= nz; ++k) {
    ZO(z0, k)
    for (i = 1; i <= nx; ++i) {
      XO(x0, i)
      exact_(i, 1, k, &UA[5 * (x0 + y1 + z0)]);
      exact_(i, ny, k, &UA[5 * (x0 + yn + z0)]);
    }
  }
  
/* ***set the dependent variable values along east and west faces */
  
  for (k = 1; k <= nz; ++k) {
    ZO(z0, k)
    for (j = 1; j <= ny; ++j) {
      YO(y0, j)
      exact_(1, j, k, &UA[5 * (x1 + y0 + z0)]);
      exact_(nx, j, k, &UA[5 * (xn + y0 + z0)]);
    }
  }
  return(0);
} /* setbv_ */


int setiv_()
{
  double peta, zeta;
  int i, j, k, m;
  double pzeta, xi, eta, pxi;
  int x0, x1, xn, y0, y1, yn, z0, z1, zn;
  
/* ***set the initial values of independent variables based on tri-linear */
/*   interpolation of boundary values in the computational space. */

  XO(x1, 1)
  YO(y1, 1)
  ZO(z1, 1)
  XO(xn, nx)
  YO(yn, ny)
  ZO(zn, nz)
  
  for (k = 2; k <= (nz-1); ++k) {
    ZO(z0, k)  
    zeta = (double) (k - 1) / (nz - 1);
    for (j = 2; j <= (ny-1); ++j) {
      YO(y0, j)
      eta = (double) (j - 1) / (ny - 1);
      for (i = 2; i <= (nx-1); ++i) {
	xi = (double) (i - 1) / (nx - 1);
        XO(x0, i)
	for (m = 1; m <= 5; ++m) {
	  pxi = (1. - xi) * UA[m + 5 * (x1 + y0 + z0) - 1] + xi * UA[m + 5 * (xn + y0 + z0) - 1];	  
	  peta = (1. - eta) * UA[m + 5 * (x0 + y1 + z0) - 1] + eta * UA[m + 5 * (x0 + yn + z0) - 1];
	  
	  pzeta = (1. - zeta) * UA[m + 5 * (x0 + y0 + z1) - 1] + zeta * UA[m + 5 * (x0 + y0 + zn) - 1];
	  
	  UA[m + 5 * (x0 + y0 + z0) - 1] = pxi + 
	    peta + pzeta - pxi * peta - peta * pzeta - pzeta *
	      pxi + pxi * peta * pzeta;
	}
      }
    }
  }

  return(0);
} /* setiv_ */

int badi_(dt, tolrsd, tx, ty, tz)
double dt, *tolrsd, *tx, *ty, *tz;
{
  extern int jacx_(), jacy_(), jacz_();
  extern int btridx_(), btridy_(), btridz_();
  extern int rhs_(), maxnorm_();

  clock_t tstart, tend;
  int imax[5], jmax[5], kmax[5]; 
  int idmax[5], jdmax[5], kdmax[5];
  int lnorm, istep;
  extern int l2norm_();
  double delunm[5];

  int i, j, k, m, x, y, z, offset;

/* ***to perform pseudo-time stepping iterations */
/*   for five coupled, nonlinear pde's. */
  
  lnorm = 2;
  
/* ***compute the steady-state residuals */
  
  rhs_(nx, ny, nz, tx, ty, tz, dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, dy5, dz1, dz2, dz3, dz4, dz5);

  if (lnorm == 1) {
    
    maxnorm_(&nx, &ny, &nz, imax, jmax, kmax, rsdnm);
    
    if (ipr == 1) {
      fprintf(out_file, "          Initial Residual norms\n\n");
      fprintf(out_file, "\n max-norm of steady-state residual for first pde  = %12.5E\n", rsdnm[0]);
      fprintf(out_file, "                                                 (%4d,%4d,%4d)\n", imax[0], jmax[0], kmax[0]); 
      fprintf(out_file, " max-norm of steady-state residual for second pde = %12.5E\n", rsdnm[1]);
      fprintf(out_file, "                                                 (%4d,%4d,%4d)\n", imax[1], jmax[1], kmax[1]); 
      fprintf(out_file, " max-norm of steady-state residual for third pde  = %12.5E\n", rsdnm[2]);
      fprintf(out_file, "                                                 (%4d,%4d,%4d)\n", imax[2], jmax[2], kmax[2]); 
      fprintf(out_file, " max-norm of steady-state residual for fourth pde = %12.5E\n", rsdnm[3]);
      fprintf(out_file, "                                                 (%4d,%4d,%4d)\n", imax[3], jmax[3], kmax[3]); 
      fprintf(out_file, " max-norm of steady-state residual for fifth pde  = %12.5E\n", rsdnm[4]);
      fprintf(out_file, "                                                 (%4d,%4d,%4d)\n", imax[4], jmax[4], kmax[4]); 
    }
  }
  else if (lnorm == 2) {
    l2norm_(&nx, &ny, &nz, rsdnm);
    
    if (ipr == 1) {
      fprintf(out_file, "          Initial Residual norms\n\n");
      fprintf(out_file, "\n RMS-norm of steady-state residual for first pde  = %12.5E\n", rsdnm[0]);
      fprintf(out_file, " RMS-norm of steady-state residual for second pde = %12.5E\n", rsdnm[1]);
      fprintf(out_file, " RMS-norm of steady-state residual for third pde  = %12.5E\n", rsdnm[2]);
      fprintf(out_file, " RMS-norm of steady-state residual for fourth pde = %12.5E\n", rsdnm[3]);
      fprintf(out_file, " RMS-norm of steady-state residual for fifth pde  = %12.5E\n", rsdnm[4]);
    }
  }
  
/* ***begin pseudo-time stepping iterations */
  
  tstart = clock();
  
  for (istep = 1; istep <= itmax; ++istep) {
    
    if (istep % inorm == 0 && ipr == 1)
      fprintf(out_file, "\n     pseudo-time Block ADI iteration no.=%4d\n\n", istep);
    
    for (k = 1; k <= nz; ++k) {
      ZO(z, k)
      for (j = 1; j <= ny; ++j) {
	YO(y, j)
	for (i = 1; i <= nx; ++i) {
	  XO(x, i)
	  offset = 5 * (x + y + z) - 1;
	  for (m = 1; m <= 5; ++m)
	    RSD[m + offset] = dt * RSD[m + offset];
	}
      }
    }
    
/* ***perform 3-factor, Block ADI iterations */

/* ***perform the xi-direction sweep */
    
    jacx_(dt, tx[0], tx[1], dx1, dx2, dx3, dx4, dx5, nx, ny, nz);
    btridx_(nx, ny, nz);

/* ***perform the eta-direction sweep */
    
    jacy_(dt, ty[0], ty[1], dy1, dy2, dy3, dy4, dy5, nx, ny, nz);
    btridy_(nx, ny, nz);
    
/* ***perform the zeta-direction sweep */

    jacz_(dt, tz[0], tz[1], dz1, dz2, dz3, dz4, dz5, nx, ny, nz);
    btridz_(nx, ny, nz);

/* ***update the variables */
    
    for (k = 1; k <= nz; ++k) {
      ZO(z, k)
      for (j = 1; j <= ny; ++j) {
	YO(y, j)
	for (i = 1; i <= nz; ++i) {
	  XO(x, i)
	  offset = 5 * (x + y + z) - 1;
	  for (m = 1; m <= 5; ++m)
	    UA[m + offset] += RSD[m + offset];
	}
      }
    }

    /* ***compute the max-norms of pseudo-time iteration corrections */
    
    if (istep % inorm == 0) {
      
      if (lnorm == 1) {
	
	maxnorm_(&nx, &ny, &nz, idmax, jdmax, kdmax, delunm);
	
	if (ipr == 1) {
	  fprintf(out_file, "\n max-norm of Block ADI-iteration correction for first pde  = %12.5E\n", delunm[0]);
	  fprintf(out_file, "                                                          (%4d,%4d,%4d)", idmax[0], jdmax[0], kdmax[0]);
	  fprintf(out_file, " max-norm of Block ADI-iteration correction for second pde = %12.5E\n", delunm[1]);
	  fprintf(out_file, "                                                          (%4d,%4d,%4d)", idmax[1], jdmax[1], kdmax[1]);
	  fprintf(out_file, " max-norm of Block ADI-iteration correction for third pde  = %12.5E\n", delunm[2]);
	  fprintf(out_file, "                                                          (%4d,%4d,%4d)", idmax[2], jdmax[2], kdmax[2]);
	  fprintf(out_file, " max-norm of Block ADI-iteration correction for fourth pde = %12.5E\n", delunm[3]);
	  fprintf(out_file, "                                                          (%4d,%4d,%4d)", idmax[3], jdmax[3], kdmax[3]);
	  fprintf(out_file, " max-norm of Block ADI-iteration correction for fifth pde  = %12.5E\n", delunm[4]);
	  fprintf(out_file, "                                                          (%4d,%4d,%4d)", idmax[4], jdmax[4], kdmax[4]);
	}
      }
      else if (lnorm == 2) {
	
	l2norm_(&nx, &ny, &nz, delunm);
	
	if (ipr == 1) {
	  fprintf(out_file, "\n RMS-norm of Block ADI-iteration correction for first pde  = %12.5E\n", delunm[0]);
	  fprintf(out_file, " RMS-norm of Block ADI-iteration correction for second pde = %12.5E\n", delunm[1]);
	  fprintf(out_file, " RMS-norm of Block ADI-iteration correction for third pde  = %12.5E\n", delunm[2]);
	  fprintf(out_file, " RMS-norm of Block ADI-iteration correction for fourth pde = %12.5E\n", delunm[3]);
	  fprintf(out_file, " RMS-norm of Block ADI-iteration correction for fifth pde  = %12.5E\n", delunm[4]);
	}
      }
    }
    
/* ***compute the steady-state residuals */
    
  rhs_(nx, ny, nz, tx, ty, tz, dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, dy5, dz1, dz2, dz3, dz4, dz5);
    
    if (istep % inorm == 0 || istep == itmax) {
      
      if (lnorm == 1) {
	
	maxnorm_(&nx, &ny, &nz, imax, jmax, kmax, rsdnm);
	
	if (ipr == 1) {
	  fprintf(out_file, "\n max-norm of steady-state residual for first pde  = %12.5E\n", rsdnm[0]);
	  fprintf(out_file, "                                                 (%4d,%4d,%4d)\n", imax[0], jmax[0], kmax[0]); 
	  fprintf(out_file, " max-norm of steady-state residual for second pde = %12.5E\n", rsdnm[1]);
	  fprintf(out_file, "                                                 (%4d,%4d,%4d)\n", imax[1], jmax[1], kmax[1]); 
	  fprintf(out_file, " max-norm of steady-state residual for third pde  = %12.5E\n", rsdnm[2]);
	  fprintf(out_file, "                                                 (%4d,%4d,%4d)\n", imax[2], jmax[2], kmax[2]); 
	  fprintf(out_file, " max-norm of steady-state residual for fourth pde = %12.5E\n", rsdnm[3]);
	  fprintf(out_file, "                                                 (%4d,%4d,%4d)\n", imax[3], jmax[3], kmax[3]); 
	  fprintf(out_file, " max-norm of steady-state residual for fifth pde  = %12.5E\n", rsdnm[4]);
	  fprintf(out_file, "                                                 (%4d,%4d,%4d)\n", imax[4], jmax[4], kmax[4]); 
	}
      } 
      else if (lnorm == 2) {
	
	l2norm_(&nx, &ny, &nz, rsdnm);
	
	if (ipr == 1) {
	  fprintf(out_file, "\n RMS-norm of steady-state residual for first pde  = %12.5E\n", rsdnm[0]);
	  fprintf(out_file, " RMS-norm of steady-state residual for second pde = %12.5E\n", rsdnm[1]);
	  fprintf(out_file, " RMS-norm of steady-state residual for third pde  = %12.5E\n", rsdnm[2]);
	  fprintf(out_file, " RMS-norm of steady-state residual for fourth pde = %12.5E\n", rsdnm[3]);
	  fprintf(out_file, " RMS-norm of steady-state residual for fifth pde  = %12.5E\n", rsdnm[4]);
	}
      }
    }
    
/****check the pseudo-time iteration residuals against the tolerance levels*/
    
    if (rsdnm[0] < tolrsd[0] && rsdnm[1] < 
	tolrsd[1] && rsdnm[2] < tolrsd[2] 
	&& rsdnm[3] < tolrsd[3] && rsdnm[4]
	< tolrsd[4])
      
      fprintf(out_file, "\n convergence was achieved after %4d pseudo-time steps\n", istep);
    
  }
  
  tend = clock();
  
  looptotal = tend - tstart;

  looptotal = 0;
  
  return(0);
  
} /* badi_ */


int maxnorm_(nx, ny, nz, imax, jmax, kmax, vnm)
int *nx, *ny, *nz, *imax, *jmax, *kmax;
double *vnm;
{
  int i, j, k, m, x, y, z, offset;
  double t1;
  
/* ***to compute the max-norm of vector v. */
  
/* Parameter adjustments */
  --vnm;
  --kmax;
  --jmax;
  --imax;
  
  for (m = 1; m <= 5; ++m)
    vnm[m] = -1e10;
  
  for (k = 1; k <= *nz; ++k) {
    ZO(z, k)
    for (j = 1; j <= *ny; ++j) {
      YO(y, j)
      for (i = 1; i <= *nx; ++i) {
	XO(x, i)
        offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m) {
	  t1 = abs(RSD[m + offset]);
	  if (vnm[m] < t1) {
	    vnm[m] = t1;
	    imax[m] = i;
	    jmax[m] = j;
	    kmax[m] = k;
	  }
	}
      }
    }
  }
  return(0);
} /* maxnorm_ */


int l2norm_(nx, ny, nz, sum)
int *nx, *ny, *nz;
double *sum;
{
  double sqrt();
  int i, j, k, m, x, y, z, offset;
  
/* ***to compute the l2-norm of vector v. */
  
/* Parameter adjustments */
  --sum;

  for (m = 1; m <= 5; ++m)
    sum[m] = 0.0;
  
  for (k = 2; k <= *nz - 1; ++k) {
    ZO(z, k)
    for (j = 2; j <= *ny - 1; ++j) {
      YO(y, j)
      for (i = 2; i <= *nx - 1; ++i) {
	XO(x, i)
        offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m)
	  sum[m] += RSD[m + offset] * RSD[m + offset];
      }
    }
  }
  
  for (m = 1; m <= 5; ++m)
    sum[m] = sqrt(sum[m] / ((*nx - 2) * (*ny - 2) * (*nz - 2)));
  
  return(0);
} /* l2norm_ */



int exact_(i, j, k, u000ijk)
int i, j, k;
double *u000ijk;
{
  double zeta;
  int m;
  double xi, eta;
  
/* ***compute the exact solution at (i,j,k) */
  
/* Function Body */
  xi = (double) (i - 1) / (nx - 1);
  eta = (double) (j - 1) / (ny - 1);
  zeta = (double) (k - 1) / (nz - 1);
  
  for (m = 1; m <= 5; ++m) {
    
    u000ijk[m - 1] = CE[m - 1] + CE[m + 4] * xi + 
               CE[m + 9] * eta + CE[m + 14] * zeta + 
               CE[m + 19] * xi * xi + CE[m + 24] * eta * 
               eta + CE[m + 29] * zeta * zeta + CE[m + 34] 
               * xi * xi * xi + CE[m + 39] * eta * eta * eta + 
               CE[m + 44] * zeta * zeta * zeta + CE[m + 49]
               * xi * xi * xi * xi + CE[m + 54] * eta * eta * eta *
               eta + CE[m + 59] * zeta * zeta * zeta * zeta;
  }
  return(0);
} /* exact_ */

int error_(errnm)
double *errnm;
{
  
/* System generated locals */
  double d_1;
  
/* Builtin functions */
  double sqrt();
  
/* Local variables */
  int imax[5], jmax[5], kmax[5], i, j, k, m;
  extern int exact_();
  int lnorm;
  double u000ijk[5], errmax[5], tmp;
  int x, y, z;
  
/* ***compute the solution error */
  
  lnorm = 2;
  
  if (lnorm == 1) {
    
    for (m = 1; m <= 5; ++m) {
      errmax[m - 1] = -1e20;
    }
    
    for (k = 2; k <= nz - 1; ++k) {
      ZO(z, k)  
      for (j = 2; j <= ny - 1; ++j) {
        YO(y, j)
	for (i = 2; i <= nx - 1; ++i) {
          XO(x, i)
	  
	  exact_(i, j, k, u000ijk);
	  
	  for (m = 1; m <= 5; ++m) {
	    
	    tmp = (d_1 = u000ijk[m - 1] - UA[m + 5 * (x + y + z) - 1], abs(d_1));
	    
	    if (tmp > errmax[m - 1]) {
	      
	      errmax[m - 1] = tmp;
	      imax[m - 1] = i;
	      jmax[m - 1] = j;
	      kmax[m - 1] = k;
	    }
	  }
	}
      }
    }
    
    fprintf(out_file, "     max. error in soln. to first pde  = %12.4E\n",errmax[0]);
    fprintf(out_file, "     and its location                  = (%4d,%4d,%4d)\n", imax[0], jmax[0], kmax[0]);
    fprintf(out_file, "     max. error in soln. to second pde = %12.4E\n",errmax[1]);
    fprintf(out_file, "     and its location                  = (%4d,%4d,%4d)\n", imax[1], jmax[1], kmax[1]);
    fprintf(out_file, "     max. error in soln. to third pde  = %12.4E\n",errmax[2]);
    fprintf(out_file, "     and its location                  = (%4d,%4d,%4d)\n", imax[2], jmax[2], kmax[2]);
    fprintf(out_file, "     max. error in soln. to fourth pde = %12.4E\n",errmax[3]);
    fprintf(out_file, "     and its location                  = (%4d,%4d,%4d)\n", imax[3], jmax[3], kmax[3]);
    fprintf(out_file, "     max. error in soln. to fifth pde  = %12.4E\n",errmax[4]);
    fprintf(out_file, "     and its location                  = (%4d,%4d,%4d)\n", imax[4], jmax[4], kmax[4]);
    
  } 
  else if (lnorm == 2) {
    
    for (m = 1; m <= 5; ++m) {
      errnm[m - 1] = 0.;
    }
    
    for (k = 2; k <= nz - 1; ++k) {
      ZO(z, k)  
      for (j = 2; j <= ny - 1; ++j) {
        YO(y, j)
	for (i = 2; i <= nx - 1; ++i) {
          XO(x, i)
	  
	  exact_(i, j, k, u000ijk);
	  for (m = 1; m <= 5; ++m) {
	    
	    tmp = u000ijk[m - 1] - UA[m + 5 * (x + y + z) - 1];
	    
	    /* Computing 2nd power */
	    d_1 = tmp;
	    errnm[m - 1] += d_1 * d_1;
	  }
	}
      }
    }
    
    for (m = 1; m <= 5; ++m) {
      errnm[m - 1] = sqrt(errnm[m - 1] / ((nx - 2) * (ny - 2) * (nz - 2)));
    }
    
    fprintf(out_file, "\n RMS-norm of error in soln. to first pde  = %12.5E\n", errnm[0]);
    fprintf(out_file, " RMS-norm of error in soln. to second pde = %12.5E\n", errnm[1]);
    fprintf(out_file, " RMS-norm of error in soln. to third pde  = %12.5E\n", errnm[2]);
    fprintf(out_file, " RMS-norm of error in soln. to fourth pde = %12.5E\n", errnm[3]);
    fprintf(out_file, " RMS-norm of error in soln. to fifth pde  = %12.5E\n", errnm[4]);
  }
  return(0);
} /* error_ */

int pintgr_(frc, dxi, deta, dzeta)
double *frc, dxi, deta, dzeta;
{
  
/* System generated locals */
  double d_1, d_2, d_3;
  
/* Local variables */
  int i, j, k;
  double frc1, frc2, frc3, phi1[S1*S2], phi2[S1*S2];
  int x0, x1, xn, y0, y1, yn, z0, z1, zn;

/* ***compute the surface integral */
  
  XO(x1, 2)
  YO(y1, 2)
  ZO(z1, 3)
  XO(xn, (nx-1))
  YO(yn, (ny-2))
  ZO(zn, (nz-1))

  for (j = 2; j <= ny - 2; ++j) {
    YO(y0, j)
    for (i = 2; i <= nx - 1; ++i) {
      XO(x0, i)
/* Computing 2nd power */
      d_1 = UA[2 + 5 * (x0 + y0 + z1) - 1];
/* Computing 2nd power */
      d_2 = UA[3 + 5 * (x0 + y0 + z1) - 1];
/* Computing 2nd power */
      d_3 = UA[4 + 5 * (x0 + y0 + z1) - 1];
      phi1[i + j * S1 - (S1+1)] = (UA[5 + 5 * (x0 + y0 + z1) - 1] - (d_1 * d_1 + d_2 * d_2 + d_3 * d_3) * .5 
                    / UA[1 + 5 * (x0 + y0 + z1) - 1]) *
                    .4;
      
/* Computing 2nd power */
      d_1 = UA[2 + 5 * (x0 + y0 + zn) - 1];
/* Computing 2nd power */
      d_2 = UA[3 + 5 * (x0 + y0 + zn) - 1];
/* Computing 2nd power */
      d_3 = UA[4 + 5 * (x0 + y0 + zn) - 1];
      phi2[i + j * S1 - (S1+1)] = (UA[5 + 5 * (x0 + y0 + zn) - 1] - (d_1 * d_1 + d_2 * d_2 + d_3 * d_3) * .5 
                    / UA[1 + 5 * (x0 + y0 + zn) - 1]) *
                    .4;
    }
  }
  frc1 = 0.;
  
  for (j = 2; j <= ny - 3; ++j) {
    for (i = 2; i <= nx - 2; ++i) {
      
      frc1 += phi1[i + j * S1 - (S1+1)] + phi1[i + 1 + j * S1 - (S1+1)] + phi1[
                    i + (j + 1) * S1 - (S1+1)] + phi1[i + 1 + (j + 1) * S1 - (S1+1)] 
                    + phi2[i + j * S1 - (S1+1)] + phi2[i + 1 + j * S1 - (S1+1)] + 
                    phi2[i + (j + 1) * S1 - (S1+1)] + phi2[i + 1 + (j + 1) * S1 - 
                    (S1+1)];
    }
  }
  
  frc1 = dxi * deta * frc1;
  
  for (k = 3; k <= nz - 1; ++k) {
    ZO(z0, k)
    for (i = 2; i <= nx - 1; ++i) {
      XO(x0, i)
/* Computing 2nd power */
      d_1 = UA[2 + 5 * (x0 + y1 + z0) - 1];
/* Computing 2nd power */
      d_2 = UA[3 + 5 * (x0 + y1 + z0) - 1];
/* Computing 2nd power */
      d_3 = UA[4 + 5 * (x0 + y1 + z0) - 1];
      phi1[i + k * S1 - (S1+1)] = (UA[5 + 5 * (x0 + y1 + z0) - 1] - (d_1 * d_1 + d_2 * d_2 + d_3 * d_3) * .5 
                    / UA[1 + 5 * (x0 + y1 + z0) - 1]) *
                    .4;
      
/* Computing 2nd power */
      d_1 = UA[2 + 5 * (x0 + yn + z0) - 1];
/* Computing 2nd power */
      d_2 = UA[3 + 5 * (x0 + yn + z0) - 1];
/* Computing 2nd power */
      d_3 = UA[4 + 5 * (x0 + yn + z0) - 1];
      phi2[i + k * S1 - (S1+1)] = (UA[5 + 5 * (x0 + yn + z0) - 1] - (d_1 * d_1 + d_2 * d_2 + d_3 * d_3) * .5 
                    / UA[1 + 5 * (x0 + yn + z0) - 1]) *
                    .4;
    }
  }
  
  frc2 = 0.;
  
  for (k = 3; k <= nz - 2; ++k) {
    for (i = 2; i <= nx - 2; ++i) {
      frc2 += phi1[i + k * S1 - (S1+1)] + phi1[i + 1 + k * S1 - (S1+1)] + phi1[
                    i + (k + 1) * S1 - (S1+1)] + phi1[i + 1 + (k + 1) * S1 - (S1+1)] 
                    + phi2[i + k * S1 - (S1+1)] + phi2[i + 1 + k * S1 - (S1+1)] + 
                    phi2[i + (k + 1) * S1 - (S1+1)] + phi2[i + 1 + (k + 1) * S1 - 
                    (S1+1)];
    }
  }
  
  frc2 = dxi * dzeta * frc2;
  
  for (k = 3; k <= nz - 1; ++k) {
    ZO(z0, k)
    for (j = 2; j <= ny - 2; ++j) {
      YO(y0, j)
/* Computing 2nd power */
      d_1 = UA[2 + 5 * (x1 + y0 + z0) - 1];
/* Computing 2nd power */
      d_2 = UA[3 + 5 * (x1 + y0 + z0) - 1];
/* Computing 2nd power */
      d_3 = UA[4 + 5 * (x1 + y0 + z0) - 1];
      phi1[j + k * S1 - (S1+1)] = (UA[5 + 5 * (x1 + y0 + z0) - 1] - (d_1 * d_1 + d_2 * d_2 + d_3 * d_3) * .5 
                    / UA[1 + 5 * (x1 + y0 + z0) - 1]) *
                    .4;
      
/* Computing 2nd power */
      d_1 = UA[2 + 5 * (xn + y0 + z0) - 1];
/* Computing 2nd power */
      d_2 = UA[3 + 5 * (xn + y0 + z0) - 1];
/* Computing 2nd power */
      d_3 = UA[4 + 5 * (xn + y0 + z0) - 1];
      phi2[j + k * S1 - (S1+1)] = (UA[5 + 5 * (xn + y0 + z0) - 1] - (d_1 * d_1 + d_2 * d_2 + d_3 * d_3) * .5 
                    / UA[1 + 5 * (xn + y0 + z0) - 1]) *
                    .4;
    }
  }
  frc3 = 0.;

  for (k = 3; k <= nz - 2; ++k) {
    for (j = 2; j <= ny - 3; ++j) {
      
      frc3 += phi1[j + k * S1 - (S1+1)] + phi1[j + 1 + k * S1 - (S1+1)] + phi1[
                    j + (k + 1) * S1 - (S1+1)] + phi1[j + 1 + (k + 1) * S1 - (S1+1)] 
                    + phi2[j + k * S1 - (S1+1)] + phi2[j + 1 + k * S1 - (S1+1)] + 
                    phi2[j + (k + 1) * S1 - (S1+1)] + phi2[j + 1 + (k + 1) * S1 - 
                    (S1+1)];
    }
  }
  
  frc3 = deta * dzeta * frc3;
  *frc = (frc1 + frc2 + frc3) * .25;
  fprintf(out_file, "\n\n     surface integral = %12.5E\n\n", *frc);
  return(0);
} /* pintgr_ */


int verify_(xcr, xce, xci)
double *xcr, *xce, *xci;
{
  
/* System generated locals */
  double d_1;
  
/* Local variables */
  int m;
  double xre[5], tmp, xri, xrr[5], epsilon;
  
/* ***verification routine */
  
/* ***tolerance level */
  
/* Parameter adjustments */
  --xce;
  --xcr;
  
/* Function Body */
  epsilon = (double)1e-8;
  
  if (nx == 12 && ny == 12 && nz == 12) {
    
/***Reference values of RMS-norms of residual, for the (12X12X12) grid, */
/*   after 60 time steps, with  DT = 1.0e-02 */
    
    xrr[0] = .17034283709541311;
    xrr[1] = .012975252070034097;
    xrr[2] = .032527926989486055;
    xrr[3] = .026436421275166801;
    xrr[4] = .1921178413174443;
    
/***Reference values of RMS-norms of solution error, for the (12X12X12) grid,*/
/*   after 60 time steps, with  DT = 1.0e-02 */
    
    xre[0] = 4.9976913345811579e-4;
    xre[1] = 4.5195666782961927e-5;
    xre[2] = 7.3973765172921357e-5;
    xre[3] = 7.3821238632439731e-5;
    xre[4] = 8.9269630987491446e-4;
    
/* ***Reference value of surface integral, for the (12X12X12) grid, */
/*   after 60 time steps, with  DT = 1.0e-02 */
    
    xri = (double)7.8414305521809488;
    
/* ***verification test for residuals */
    
    for (m = 1; m <= 5; ++m) {
      tmp = (d_1 = (xcr[m] - xrr[m - 1]) / xrr[m - 1], abs(d_1));
      if (tmp > epsilon) {
	fprintf(out_file, "\n     VERIFICATION TEST FOR RESIDUALS FAILED\n");
	goto L100;
      }
    }
    
    fprintf(out_file, "\n     VERIFICATION TEST FOR RESIDUALS IS SUCCESSFUL\n");
    
/* ***verification test for solution error */
    
  L100:
    for (m = 1; m <= 5; ++m) {
      tmp = (d_1 = (xce[m] - xre[m - 1]) / xre[m - 1], abs(d_1));
      if (tmp > epsilon) {
	fprintf(out_file, "\n     VERIFICATION TEST FOR SOLUTION ERRORS FAILED\n");
	goto L200;
      }
    }
    
    fprintf(out_file, "\n     VERIFICATION TEST FOR SOLUTION ERRORS IS SUCCESSFUL\n");

/* ***verification test for surface integral */
    
  L200:
    
    tmp = (d_1 = (*xci - xri) / xri, abs(d_1));
    if (tmp > epsilon)
      fprintf(out_file, "\n     VERIFICATION TEST FOR SURFACE INTEGRAL FAILED\n");
    else
      fprintf(out_file, "\n     VERIFICATION TEST FOR SURFACE INTEGRAL IS SUCCESSFUL\n");
    
    fprintf(out_file, "\n          CAUTION\n");
    fprintf(out_file, "\n     REFERENCE VALUES CURRENTLY IN THIS VERIFICATION ROUTINE\n");
    fprintf(out_file, "     ARE VALID ONLY FOR RUNS WITH THE FOLLOWING PARAMETER VALUES:\n");
    fprintf(out_file, "\n     NX = 12;  NY = 12;  NZ = 12\n\n     ITMAX = 60\n");
    fprintf(out_file, "\n     DT = 1.0e-02\n");
    fprintf(out_file, "\n     CHANGE IN ANY OF THE ABOVE VALUES RENDER THE REFERENCE VALUES\n");
    fprintf(out_file, "     INVALID AND CAUSES A FAILURE OF THE VERIFICATION TEST.\n");
    
  } 
  else if (nx == 64 && ny == 64 && nz == 64) {
    
/* ***Reference values of RMS-norms of residual, for the (64X64X64) gr
       id, */
/*   after 200 time steps, with  DT = 1.0e-03 */
    
    xrr[0] = 108.06346714637264;
    xrr[1] = 11.319730901220813;
    xrr[2] = 25.974354511582465;
    xrr[3] = 23.66562254467891;
    xrr[4] = 252.78963211748344;
    
/****Reference values of RMS-norms of solution error, for the (64X64X6
      4) grid,*/
/*   after 200 time steps, with  DT = 1.0e-03 */
    
    xre[0] = 4.2348416040525025;
    xre[1] = .44390282496995698;
    xre[2] = .9669248013634565;
    xre[3] = .88302063039765474;
    xre[4] = 9.7379901770829278;
    
/* ***Reference value of surface integral, for the (64X64X64) grid, */
    
/*   after 200 time steps, with DT = 1.0e-03 */
    
    xri = 14.101360012200832;
    
/* ***verification test for residuals */
    
    for (m = 1; m <= 5; ++m) {
      
      tmp = (d_1 = (xcr[m] - xrr[m - 1]) / xrr[m - 1], abs(d_1));
      
      if (tmp > epsilon) {
	
	fprintf(out_file, "\n     VERIFICATION TEST FOR RESIDUALS FAILED\n");
	goto L400;
      }
    }
    fprintf(out_file, "\n     VERIFICATION TEST FOR RESIDUALS IS SUCCESSFUL\n");
    
/* ***verification test for solution error */
    
  L400:
    
    for (m = 1; m <= 5; ++m) {
      tmp = (d_1 = (xce[m] - xre[m - 1]) / xre[m - 1], abs(d_1));
      if (tmp > epsilon) {
	
	fprintf(out_file, "\n     VERIFICATION TEST FOR SOLUTION ERRORS FAILED\n");
	goto L500;
      }
    }
    
    fprintf(out_file, "\n     VERIFICATION TEST FOR SOLUTION ERRORS IS SUCCESSFUL\n");
    
/* ***verification test for surface integral */
    
  L500:
    
    tmp = (d_1 = (*xci - xri) / xri, abs(d_1));
    
    if (tmp > epsilon) 
      fprintf(out_file, "\n     VERIFICATION TEST FOR SURFACE INTEGRAL FAILED\n");
    
    else
      fprintf(out_file, "\n     VERIFICATION TEST FOR SURFACE INTEGRAL IS SUCCESSFUL\n");
    
    fprintf(out_file, "\n          CAUTION\n");
    fprintf(out_file, "\n     REFERENCE VALUES CURRENTLY IN THIS VERIFICATION ROUTINE\n");
    fprintf(out_file, "     ARE VALID ONLY FOR RUNS WITH THE FOLLOWING PARAMETER VALUES:\n");
    fprintf(out_file, "\n     NX = 64;  NY = 64;  NZ = 64\n\n     ITMAX = 200\n");
    fprintf(out_file, "\n     DT = 8.0e-04\n");
    fprintf(out_file, "\n     CHANGE IN ANY OF THE ABOVE VALUES RENDER THE REFERENCE VALUES\n");
    fprintf(out_file, "     INVALID AND CAUSES A FAILURE OF THE VERIFICATION TEST.\n");
  } 
  else {
    fprintf(out_file, "\n FOR THE PROBLEM PARAMETERS IN USE NO REFERENCE VALUES ARE PROVIDED\n");
    fprintf(out_file, " IN THE CURRENT VERIFICATION ROUTINE - NO VERIFICATION TEST WAS PERFORMED\n");
  }
  
  return(0);
} /* verify_ */

/* Subroutine for cubic decomposition, breaks domain into n cubes */
void partition_domain(domain, cube, num_procs, my_proc)
int *domain, *cube, num_procs, my_proc;
{
  int counter=0;

  *(domain) = 1;
  *(domain + 1) = 1;
  *(domain + 2) = 1;

  while (num_procs >= 2) {
    *(domain + counter) *= 2;
    counter = (counter + 1) % 3;
    num_procs /= 2;
  }
}

main()
{
  char buf[MAXLINE];
  extern int erhs_(), setbv_(), setiv_(), error_(), 
  pintgr_(), verify_(), badi_();
  clock_t tstart, tend;
  int x_cube, y_cube, z_cube;
  int temp1, temp2, temp3, temp4;
  double frc, dt, tolrsd[5], errnm[5]; 
  double dxi, deta, dzeta, tx[3], ty[3], tz[3];
  char *prog_name = "appbt";
  char *in_name = "appbt.input";
  char *out_name = "output.data";
  int iout, i;

/* ***driver for the performance evaluation of the solver for */
/*   five coupled, nonlinear partial differential equations. */
  
/* Open file for input according to predefined name */
  if ((in_file = fopen(in_name, "r")) == NULL) {
    fprintf(stderr, "%s: can't open %s\n", prog_name, in_name);
    exit(1);
  }
  
/* ***read the unit number for output data */
  
  fgets(buf, MAXLINE, in_file);  /* Read blank lines in between input */
/*  fgets(buf, MAXLINE, in_file);  */ /* The input file doesn't have this */
  fgets(buf, MAXLINE, in_file);       /* blank line but the fortran version */
  sscanf(buf, "%d\n", &iout);         /* does.  So it's commented out */
  
/* ***flag that controls printing of the progress of iterations */
  
  fgets(buf, MAXLINE, in_file);  
  fgets(buf, MAXLINE, in_file);  
  fgets(buf, MAXLINE, in_file);  
  sscanf(buf, "%d%d\n", &ipr, &inorm);
  
/* ***set the maximum number of pseudo-time steps to be taken */
  
  fgets(buf, MAXLINE, in_file);  
  fgets(buf, MAXLINE, in_file);  
  fgets(buf, MAXLINE, in_file);  
  sscanf(buf, "%d\n", &itmax);
  
/* ***set the magnitude of the time step */
  
  fgets(buf, MAXLINE, in_file);  
  fgets(buf, MAXLINE, in_file);  
  fgets(buf, MAXLINE, in_file);  
  sscanf(buf, "%lf\n", &dt);
  
/* ***set the value of over-relaxation factor for SSOR iterations */
  
  fgets(buf, MAXLINE, in_file);  
  fgets(buf, MAXLINE, in_file);  
  fgets(buf, MAXLINE, in_file);  
  sscanf(buf, "%lf\n", &tolrsd[0]); /* This was omega, but the value is unused */
  
/* ***set the steady-state residual tolerance levels */
  
  fgets(buf, MAXLINE, in_file);  
  fgets(buf, MAXLINE, in_file);  
  fgets(buf, MAXLINE, in_file);  
  sscanf(buf, "%lf%lf%lf%lf%lf\n", &tolrsd[0], &tolrsd[1], &tolrsd[2], &tolrsd[3], &tolrsd[4]);
  
/* ***read problem specification parameters */
/* ***specify the number of grid points in xi, eta and zeta directions */  
  
  fgets(buf, MAXLINE, in_file);  
  fgets(buf, MAXLINE, in_file);  
  fgets(buf, MAXLINE, in_file);
  fclose(in_file);
  sscanf(buf, "%d%d%d\n", &nx, &ny, &nz);
  
/* ***open the file for output data */
  
  if (iout == 7) {
    if ((out_file = fopen(out_name, "w")) == NULL)
      fprintf(stderr, "%s: can't open %s\n", prog_name, out_name);
    exit(1);
  }
  else
    out_file = stdout;

  if (nx < 5 || ny < 5 || nz < 5) {
    fprintf(stderr, "PROBLEM SIZE IS TOO SMALL -\nSET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5\n");
    exit(1);
  }
  
  if (nx > S1 || ny > S2 || nz > S3) {
    fprintf(stderr, "PROBLEM SIZE IS TOO LARGE -\nNX, NY AND NZ SHOULD BE LESS THAN OR EQUAL TO ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY\n");
    exit(1);
  }
  
  tstart = clock();

  dxi = 1. / (nx - 1);
  deta = 1. / (ny - 1);
  dzeta = 1. / (nz - 1);
  
  tx[0] = 1. / (dxi * dxi);
  tx[1] = 1. / (dxi * 2.);
  tx[2] = 1. / dxi;
  
  ty[0] = 1. / (deta * deta);
  ty[1] = 1. / (deta * 2.);
  ty[2] = 1. / deta;
  
  tz[0] = 1. / (dzeta * dzeta);
  tz[1] = 1. / (dzeta * 2.);
  tz[2] = 1. / dzeta;
  
/* Set the diffusion constants */

  dx1 = .75;
  dx2 = dx1;
  dx3 = dx1;
  dx4 = dx1;
  dx5 = dx1;
  
  dy1 = .75;
  dy2 = dy1;
  dy3 = dy1;
  dy4 = dy1;
  dy5 = dy1;
  
  dz1 = 1.;
  dz2 = dz1;
  dz3 = dz1;
  dz4 = dz1;
  dz5 = dz1;
  
  GMEM = (Global_Struct *) malloc(sizeof(Global_Struct));

/* ***coefficients of the exact solution to the first pde */
  
  CE[0] = 2.;
  CE[5] = 0.;
  CE[10] = 0.;
  CE[15] = 4.;
  CE[20] = 5.;
  CE[25] = 3.;
  CE[30] = .5;
  CE[35] = .02;
  CE[40] = .01;
  CE[45] = .03;
  CE[50] = .5;
  CE[55] = .4;
  CE[60] = .3;
  
/* ***coefficients of the exact solution to the second pde */
  
  CE[1] = 1.;
  CE[6] = 0.;
  CE[11] = 0.;
  CE[16] = 0.;
  CE[21] = 1.;
  CE[26] = 2.;
  CE[31] = 3.;
  CE[36] = .01;
  CE[41] = .03;
  CE[46] = .02;
  CE[51] = .4;
  CE[56] = .3;
  CE[61] = .5;
  
/* ***coefficients of the exact solution to the third pde */
  
  CE[2] = 2.;
  CE[7] = 2.;
  CE[12] = 0.;
  CE[17] = 0.;
  CE[22] = 0.;
  CE[27] = 2.;
  CE[32] = 3.;
  CE[37] = .04;
  CE[42] = .03;
  CE[47] = .05;
  CE[52] = .3;
  CE[57] = .5;
  CE[62] = .4;
  
/* ***coefficients of the exact solution to the fourth pde */
  
  CE[3] = 2.;
  CE[8] = 2.;
  CE[13] = 0.;
  CE[18] = 0.;
  CE[23] = 0.;
  CE[28] = 2.;
  CE[33] = 3.;
  CE[38] = .03;
  CE[43] = .05;
  CE[48] = .04;
  CE[53] = .2;
  CE[58] = .1;
  CE[63] = .3;
  
/* ***coefficients of the exact solution to the fifth pde */
  
  CE[4] = 5.;
  CE[9] = 4.;
  CE[14] = 3.;
  CE[19] = 2.;
  CE[24] = .1;
  CE[29] = .4;
  CE[34] = .3;
  CE[39] = .05;
  CE[44] = .04;
  CE[49] = .03;
  CE[54] = .1;
  CE[59] = .3;
  CE[64] = .2;
  
  partition_domain(subblocks, subblocks, 1, 1);

  x_cube = (nx - 1)/subblocks[0] + 1;
  y_cube = (ny - 1)/subblocks[1] + 1;
  z_cube = (nz - 1)/subblocks[2] + 1;

  temp4 = x_cube * y_cube;
  temp1 = x_cube * y_cube * z_cube - x_cube;
  temp2 = nx * y_cube * z_cube - temp4;
  temp3 = nx * ny * z_cube - x_cube * y_cube * z_cube;

/* Save offsets for cubic array references in pre-defined array */

  for (i = 1; i<=nx; i++)
    X_A[i-1] = i - 1 + ((i-1)/x_cube) * temp1;
  for (i = 1; i<=ny; i++)
    Y_A[i-1] = ((i-1)/y_cube) * temp2 + (i - 1) * x_cube;
  for (i = 1; i<=nz; i++)
    Z_A[i-1] = ((i-1)/z_cube) * temp3 + (i - 1) * temp4;

  setbv_(); /* ***set the boundary values for dependent variables */
  setiv_(); /* ***set the initial values for dependent variables */
  erhs_(nx, ny, nz, tx, ty, tz, dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, dy5, dz1, dz2, dz3, dz4, dz5);
  badi_(dt, tolrsd, tx, ty, tz); /* ***perform Block approximate factorization iterations */
  error_(errnm); /* ***compute the solution error */
  pintgr_(&frc, dxi, deta, dzeta); /* ***compute the surface integral */
  verify_(rsdnm, errnm, &frc); /* ***verification test */
  
  tend = clock();
  alltotal = tend - tstart;

/* ***print the CPU time */

  fprintf(out_file, "\n\n     Internal Loop CPU time = %12.4E Sec.\n\n", (float) (looptotal/1000000));
  fprintf(out_file, "     Total Program CPU time = %12.4E Sec.\n\n", (float) (alltotal/1000000));
  
} /* main */



