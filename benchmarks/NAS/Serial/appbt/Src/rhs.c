/* --- rhs.c ---
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

#include "appbt.h"

extern Global_Struct *GMEM;

int erhs_(nx, ny, nz, tx, ty, tz, dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, dy5, dz1, dz2, dz3, dz4, dz5)
int nx, ny, nz;
double *tx, *ty, *tz, dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, dy5, dz1, dz2, dz3, dz4, dz5;
{
/* System generated locals */
  double d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8;
  
/* Local variables */
  double zeta, flux[5*S1]	/* was [5][S1] */, u21im1, u31im1, 
  u41im1, u51im1, u21jm1, u31jm1, u41jm1, u51jm1, u21km1, u31km1, 
  u41km1, u51km1;
  double q, dsspm, u21, ue[5*S1], u31, u41, xi, eta, u21i, u31i, u41i, u51i;
  double u21j, u31j, u41j, u51j, u21k, u31k, u41k, u51k, tmp;
  double tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3;

  int i, j, k, m;
  int x, y, z, offset;
  
  tx1 = tx[0];
  tx2 = tx[1];
  tx3 = tx[2];
  ty1 = ty[0];
  ty2 = ty[1];
  ty3 = ty[2];
  tz1 = tz[0];
  tz2 = tz[1];
  tz3 = tz[2];

/* ***compute the right hand side based on exact solution */
  
  dsspm = 1. / 4.;
  
  for (k = 1; k <= nz; ++k) {
    ZO(z, k);
    for (j = 1; j <= ny; ++j) {
      YO(y, j);
      for (i = 1; i <= nx; ++i) {
        XO(x, i);
	offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m) {
	  FRCT[m + offset] = 0.;
	}
      }
    }
  }
  
/* ***xi-direction flux differences */
  
  for (k = 2; k <= nz - 1; ++k) {
    ZO(z, k);
    zeta = (double) (k - 1) / (nz - 1);
    for (j = 2; j <= ny - 1; ++j) {
      YO(y, j);
      eta = (double) (j - 1) / (ny - 1);
      for (i = 1; i <= nx; ++i) {
	xi = (double) (i - 1) / (nx - 1);
	for (m = 1; m <= 5; ++m) {
	  ue[m + i * 5 - 6] = CE[m - 1] + CE[m + 
                    4] * xi + CE[m + 9] * eta + CE[
                    m + 14] * zeta + CE[m + 19] * xi * xi + 
                    CE[m + 24] * eta * eta + CE[m + 
                    29] * zeta * zeta + CE[m + 34] * xi * xi 
                    * xi + CE[m + 39] * eta * eta * eta + 
                    CE[m + 44] * zeta * zeta * zeta + 
                    CE[m + 49] * xi * xi * xi * xi + 
                    CE[m + 54] * eta * eta * eta * eta + 
                    CE[m + 59] * zeta * zeta * zeta * zeta;
	}
	
	flux[i * 5 - 5] = ue[i * 5 - 4];
	
	u21 = ue[i * 5 - 4] / ue[i * 5 - 5];
	
	q = (ue[i * 5 - 4] * ue[i * 5 - 4] + ue[i * 5 - 3] * ue[i * 5 
                    - 3] + ue[i * 5 - 2] * ue[i * 5 - 2]) * .5 / ue[i * 5 - 5];
	
	flux[i * 5 - 4] = ue[i * 5 - 4] * u21 + (ue[i * 5 - 1] - q) * .4;
	flux[i * 5 - 3] = ue[i * 5 - 3] * u21;
	flux[i * 5 - 2] = ue[i * 5 - 2] * u21;
	flux[i * 5 - 1] = (ue[i * 5 - 1] * 1.4 - q * .4) * u21;
      }
      
      for (i = 2; i <= nx - 1; ++i) {
        XO(x, i);
	offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m)
	  FRCT[m + offset] -= 
                    tx2 * (flux[m + (i + 1) * 5 - 6] - flux[m 
                    + (i - 1) * 5 - 6]);
      }
      for (i = 2; i <= nx; ++i) {
	tmp = 1. / ue[i * 5 - 5];
	
	u21i = tmp * ue[i * 5 - 4];
	u31i = tmp * ue[i * 5 - 3];
	u41i = tmp * ue[i * 5 - 2];
	u51i = tmp * ue[i * 5 - 1];
	
	tmp = 1. / ue[(i - 1) * 5 - 5];
	
	u21im1 = tmp * ue[(i - 1) * 5 - 4];
	u31im1 = tmp * ue[(i - 1) * 5 - 3];
	u41im1 = tmp * ue[(i - 1) * 5 - 2];
	u51im1 = tmp * ue[(i - 1) * 5 - 1];
	
	flux[i * 5 - 4] = tx3 * 1.3333333333333333 * (u21i - u21im1);
	flux[i * 5 - 3] = tx3 * (u31i - u31im1);
	flux[i * 5 - 2] = tx3 * (u41i - u41im1);
	/* Computing 2nd power */
	d_1 = u21i;
	/* Computing 2nd power */
	d_2 = u31i;
	/* Computing 2nd power */
	d_3 = u41i;
	/* Computing 2nd power */
	d_4 = u21im1;
	/* Computing 2nd power */
	d_5 = u31im1;
	/* Computing 2nd power */
	d_6 = u41im1;
	/* Computing 2nd power */
	d_7 = u21i;
	/* Computing 2nd power */
	d_8 = u21im1;
	flux[i * 5 - 1] = tx3 * -.47999999999999987 * (d_1 * 
                    d_1 + d_2 * d_2 + d_3 * d_3 - (d_4 * d_4 + d_5 * d_5 
                    + d_6 * d_6)) + tx3 * .16666666666666666 * (
                    d_7 * d_7 - d_8 * d_8) + tx3 * 
                    1.9599999999999997 * (u51i - u51im1);
      }
      
      for (i = 2; i <= nx - 1; ++i) {
        XO(x, i);
	offset = 5 * (x + y + z) - 1;
	FRCT[1 + offset] += dx1 *
                    tx1 * (ue[(i - 1) * 5 - 5] - ue[i * 5 - 5] * 
                    2. + ue[(i + 1) * 5 - 5]);
	
	FRCT[2 + offset] = FRCT[2 + offset] + tx3 * .1 *
                    1. * (flux[(i + 1) * 5 - 4] - flux[i * 5 - 4]) + 
                    dx2 * tx1 * (ue[(i - 1) * 5 - 4] - ue[
                    i * 5 - 4] * 2. + ue[(i + 1) * 5 - 4]);
	
	FRCT[3 + offset] = FRCT[3 + offset] + tx3 * .1 *
                    1. * (flux[(i + 1) * 5 - 3] - flux[i * 5 - 3]) + 
                    dx3 * tx1 * (ue[(i - 1) * 5 - 3] - ue[
                    i * 5 - 3] * 2. + ue[(i + 1) * 5 - 3]);
	
	FRCT[4 + offset] = FRCT[4 + offset] + tx3 * .1 *
                    1. * (flux[(i + 1) * 5 - 2] - flux[i * 5 - 2]) + 
                    dx4 * tx1 * (ue[(i - 1) * 5 - 2] - ue[
                    i * 5 - 2] * 2. + ue[(i + 1) * 5 - 2]);
	
	FRCT[5 + offset] = FRCT[5 + offset] + tx3 * .1 *
                    1. * (flux[(i + 1) * 5 - 1] - flux[i * 5 - 1]) + 
                    dx5 * tx1 * (ue[(i - 1) * 5 - 1] - ue[
                    i * 5 - 1] * 2. + ue[(i + 1) * 5 - 1]);
      }
      
/* ***Fourth-order dissipation */
      
      for (m = 1; m <= 5; ++m) {
	
	FRCT[m + 5 * (X_A[1] + y + z) - 1] -= dsspm * 
                    (ue[m + 4] * 5. - ue[m + 9] * 4. + ue[m + 14]);
	
	FRCT[m + 5 * (X_A[2] + y + z) - 1] -= dsspm * 
                    (ue[m + 4] * -4. + ue[m + 9] * 6. - ue[m + 14] * 4. + 
                    ue[m + 19]);
      }
      
      for (i = 4; i <= nx - 3; ++i) {
        XO(x, i);
	offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m)
	  FRCT[m + offset] -= 
                    dsspm * (ue[m + (i - 2) * 5 - 6] - ue[m + (i - 1) 
                    * 5 - 6] * 4. + ue[m + i * 5 - 6] * 6. - ue[m + (
                    i + 1) * 5 - 6] * 4. + ue[m + (i + 2) * 5 - 6]);
      }
      for (m = 1; m <= 5; ++m) {
	FRCT[m + 5 * (X_A[nx-3] + y + z) - 1] -= dsspm * (ue[m + (nx - 4) * 5 - 6] - 
                    ue[m + (nx - 3) * 5 - 6] * 4. + ue[m + (
                    nx - 2) * 5 - 6] * 6. - ue[m + (nx - 
                    1) * 5 - 6] * 4.);
	
	FRCT[m + 5 * (X_A[nx-2] + y + z) - 1] -= dsspm * (ue[m + (nx - 3) * 5 - 6] - 
                    ue[m + (nx - 2) * 5 - 6] * 4. + ue[m + (
                    nx - 1) * 5 - 6] * 5.);
      }
    }
  }

  /* ***eta-direction flux differences */
  
  for (k = 2; k <= nz - 1; ++k) {
    ZO(z, k);
    zeta = (double) (k - 1) / (nz - 1);
    for (i = 2; i <= nx - 1; ++i) {
      XO(x, i);
      xi = (double) (i - 1) / (nx - 1);
      for (j = 1; j <= ny; ++j) {
	eta = (double) (j - 1) / (ny - 1);
	for (m = 1; m <= 5; ++m) {
	  ue[m + j * 5 - 6] = CE[m - 1] + CE[m + 
                    4] * xi + CE[m + 9] * eta + CE[
                    m + 14] * zeta + CE[m + 19] * xi * xi + 
                    CE[m + 24] * eta * eta + CE[m + 
                    29] * zeta * zeta + CE[m + 34] * xi * xi 
                    * xi + CE[m + 39] * eta * eta * eta + 
                    CE[m + 44] * zeta * zeta * zeta + 
                    CE[m + 49] * xi * xi * xi * xi + 
                    CE[m + 54] * eta * eta * eta * eta + 
                    CE[m + 59] * zeta * zeta * zeta * zeta;
	}
	
	flux[j * 5 - 5] = ue[j * 5 - 3];
	
	u31 = ue[j * 5 - 3] / ue[j * 5 - 5];
	
	q = (ue[j * 5 - 4] * ue[j * 5 - 4] + ue[j * 5 - 3] * ue[j * 5 
                    - 3] + ue[j * 5 - 2] * ue[j * 5 - 2]) * .5 / ue[j * 5 - 5];
	
	flux[j * 5 - 4] = ue[j * 5 - 4] * u31;
	flux[j * 5 - 3] = ue[j * 5 - 3] * u31 + (ue[j * 5 - 1] - q) * .4;
	flux[j * 5 - 2] = ue[j * 5 - 2] * u31;
	flux[j * 5 - 1] = (ue[j * 5 - 1] * 1.4 - q * .4) * u31;
      }
      
      for (j = 2; j <= ny - 1; ++j) {
	YO(y, j);
	offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m)
	  FRCT[m + offset] -= 
                    ty2 * (flux[m + (j + 1) * 5 - 6] - flux[m 
                    + (j - 1) * 5 - 6]);
      }
      for (j = 2; j <= ny; ++j) {
	
	tmp = 1. / ue[j * 5 - 5];
	
	u21j = tmp * ue[j * 5 - 4];
	u31j = tmp * ue[j * 5 - 3];
	u41j = tmp * ue[j * 5 - 2];
	u51j = tmp * ue[j * 5 - 1];
	
	tmp = 1. / ue[(j - 1) * 5 - 5];
	
	u21jm1 = tmp * ue[(j - 1) * 5 - 4];
	u31jm1 = tmp * ue[(j - 1) * 5 - 3];
	u41jm1 = tmp * ue[(j - 1) * 5 - 2];
	u51jm1 = tmp * ue[(j - 1) * 5 - 1];
	
	flux[j * 5 - 4] = ty3 * (u21j - u21jm1);
	flux[j * 5 - 3] = ty3 * 1.3333333333333333 * (u31j - u31jm1);
	flux[j * 5 - 2] = ty3 * (u41j - u41jm1);
	/* Computing 2nd power */
	d_1 = u21j;
	/* Computing 2nd power */
	d_2 = u31j;
	/* Computing 2nd power */
	d_3 = u41j;
	/* Computing 2nd power */
	d_4 = u21jm1;
	/* Computing 2nd power */
	d_5 = u31jm1;
	/* Computing 2nd power */
	d_6 = u41jm1;
	/* Computing 2nd power */
	d_7 = u31j;
	/* Computing 2nd power */
	d_8 = u31jm1;
	flux[j * 5 - 1] = ty3 * -.47999999999999987 * (d_1 * 
                    d_1 + d_2 * d_2 + d_3 * d_3 - (d_4 * d_4 + d_5 * d_5 
                    + d_6 * d_6)) + ty3 * .16666666666666666 * (
                    d_7 * d_7 - d_8 * d_8) + ty3 * 
                    1.9599999999999997 * (u51j - u51jm1);
      }
      
      for (j = 2; j <= ny - 1; ++j) {
	YO(y, j);
	offset = 5 * (x + y + z) - 1;
	
	FRCT[1 + offset] += dy1 *
                    ty1 * (ue[(j - 1) * 5 - 5] - ue[j * 5 - 5] * 
                    2. + ue[(j + 1) * 5 - 5]);
	
	FRCT[2 + offset] = FRCT[2 + offset] + ty3 * .1 *
                    1. * (flux[(j + 1) * 5 - 4] - flux[j * 5 - 4]) + 
                    dy2 * ty1 * (ue[(j - 1) * 5 - 4] - ue[
                    j * 5 - 4] * 2. + ue[(j + 1) * 5 - 4]);
                    
	FRCT[3 + offset] = FRCT[3 + offset] + ty3 * .1 *
                    1. * (flux[(j + 1) * 5 - 3] - flux[j * 5 - 3]) + 
                    dy3 * ty1 * (ue[(j - 1) * 5 - 3] - ue[
                    j * 5 - 3] * 2. + ue[(j + 1) * 5 - 3]);
	
	FRCT[4 + offset] = FRCT[4 + offset] + ty3 * .1 *
                    1. * (flux[(j + 1) * 5 - 2] - flux[j * 5 - 2]) + 
                    dy4 * ty1 * (ue[(j - 1) * 5 - 2] - ue[
                    j * 5 - 2] * 2. + ue[(j + 1) * 5 - 2]);
	
	FRCT[5 + offset] = FRCT[5 + offset] + ty3 * .1 *
                    1. * (flux[(j + 1) * 5 - 1] - flux[j * 5 - 1]) + 
                    dy5 * ty1 * (ue[(j - 1) * 5 - 1] - ue[
                    j * 5 - 1] * 2. + ue[(j + 1) * 5 - 1]);
      }
      
/* ***fourth-order dissipation */
      
      for (m = 1; m <= 5; ++m) {
	
	FRCT[m + 5 * (x + Y_A[1] + z) - 1] -= dsspm * 
	  (ue[m + 4] * 5. - ue[m + 9] * 4. + ue[m + 14]);
	
	FRCT[m + 5 * (x + Y_A[2] + z) - 1] -= dsspm * 
	  (ue[m + 4] * -4. + ue[m + 9] * 6. - ue[m + 14] * 4. + 
	   ue[m + 19]);
	
      }
      
      for (j = 4; j <= ny - 3; ++j) {
	YO(y, j);
	offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m)
	  FRCT[m + offset] -= 
                    dsspm * (ue[m + (j - 2) * 5 - 6] - ue[m + (j - 1) 
                    * 5 - 6] * 4. + ue[m + j * 5 - 6] * 6. - ue[m + (
                    j + 1) * 5 - 6] * 4. + ue[m + (j + 2) * 5 - 6]);
      }
      for (m = 1; m <= 5; ++m) {
	FRCT[m + 5 * (x + Y_A[ny-3] + z) - 1] -= dsspm * (ue[m + (ny - 4) * 5 - 6] - 
                    ue[m + (ny - 3) * 5 - 6] * 4. + ue[m + (
                    ny - 2) * 5 - 6] * 6. - ue[m + (ny - 
                    1) * 5 - 6] * 4.);
	
	FRCT[m + 5 * (x + Y_A[ny-2] + z) - 1] -= dsspm * (ue[m + (ny - 3) * 5 - 6] - 
                    ue[m + (ny - 2) * 5 - 6] * 4. + ue[m + (
                    ny - 1) * 5 - 6] * 5.);
      }
    }
  }
  
/* ***zeta-direction flux differences */
  
  for (j = 2; j <= ny - 1; ++j) {
    YO(y, j);
    eta = (double) (j - 1) / (ny - 1);
    
    for (i = 2; i <= nx - 1; ++i) {
      XO(x, i);
      xi = (double) (i - 1) / (nx - 1);
      for (k = 1; k <= nz; ++k) {
	zeta = (double) (k - 1) / (nz - 1);
	
	for (m = 1; m <= 5; ++m)
	  ue[m + k * 5 - 6] = CE[m - 1] + CE[m + 
                    4] * xi + CE[m + 9] * eta + CE[
                    m + 14] * zeta + CE[m + 19] * xi * xi + 
                    CE[m + 24] * eta * eta + CE[m + 
                    29] * zeta * zeta + CE[m + 34] * xi * xi 
                    * xi + CE[m + 39] * eta * eta * eta + 
                    CE[m + 44] * zeta * zeta * zeta + 
                    CE[m + 49] * xi * xi * xi * xi + 
                    CE[m + 54] * eta * eta * eta * eta + 
                    CE[m + 59] * zeta * zeta * zeta * zeta;
	
	flux[k * 5 - 5] = ue[k * 5 - 2];
	
	u41 = ue[k * 5 - 2] / ue[k * 5 - 5];
	
	q = (ue[k * 5 - 4] * ue[k * 5 - 4] + ue[k * 5 - 3] * ue[k * 5 
                    - 3] + ue[k * 5 - 2] * ue[k * 5 - 2]) * .5 / ue[k * 5 - 5];
	
	flux[k * 5 - 4] = ue[k * 5 - 4] * u41;
	flux[k * 5 - 3] = ue[k * 5 - 3] * u41;
	flux[k * 5 - 2] = ue[k * 5 - 2] * u41 + (ue[k * 5 - 1] - q) * .4;
	flux[k * 5 - 1] = (ue[k * 5 - 1] * 1.4 - q * .4) * u41;
      }
      
      for (k = 2; k <= nz - 1; ++k) {
	ZO(z, k)
        offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m)
	  FRCT[m + offset] -= 
                    tz2 * (flux[m + (k + 1) * 5 - 6] - flux[m 
                    + (k - 1) * 5 - 6]);
      }
      for (k = 2; k <= nz; ++k) {
	
	tmp = 1. / ue[k * 5 - 5];
	
	u21k = tmp * ue[k * 5 - 4];
	u31k = tmp * ue[k * 5 - 3];
	u41k = tmp * ue[k * 5 - 2];
	u51k = tmp * ue[k * 5 - 1];
	
	tmp = 1. / ue[(k - 1) * 5 - 5];
	
	u21km1 = tmp * ue[(k - 1) * 5 - 4];
	u31km1 = tmp * ue[(k - 1) * 5 - 3];
	u41km1 = tmp * ue[(k - 1) * 5 - 2];
	u51km1 = tmp * ue[(k - 1) * 5 - 1];
	
	flux[k * 5 - 4] = tz3 * (u21k - u21km1);
	flux[k * 5 - 3] = tz3 * (u31k - u31km1);
	flux[k * 5 - 2] = tz3 * 1.3333333333333333 * (u41k - u41km1);
	/* Computing 2nd power */
	d_1 = u21k;
	/* Computing 2nd power */
	d_2 = u31k;
	/* Computing 2nd power */
	d_3 = u41k;
	/* Computing 2nd power */
	d_4 = u21km1;
	/* Computing 2nd power */
	d_5 = u31km1;
	/* Computing 2nd power */
	d_6 = u41km1;
	/* Computing 2nd power */
	d_7 = u41k;
	/* Computing 2nd power */
	d_8 = u41km1;
	flux[k * 5 - 1] = tz3 * -.47999999999999987 * (d_1 * 
                    d_1 + d_2 * d_2 + d_3 * d_3 - (d_4 * d_4 + d_5 * d_5 
                    + d_6 * d_6)) + tz3 * .16666666666666666 * (
                    d_7 * d_7 - d_8 * d_8) + tz3 * 
                    1.9599999999999997 * (u51k - u51km1);
      }
      
      for (k = 2; k <= nz - 1; ++k) {
	ZO(z, k);
	offset = 5 * (x + y + z) - 1;
	FRCT[1 + offset] += dz1 *
                    tz1 * (ue[(k + 1) * 5 - 5] - ue[k * 5 - 5] * 
                    2. + ue[(k - 1) * 5 - 5]);
	
	FRCT[2 + offset] = FRCT[2 + offset] + tz3 * .1 *
                    1. * (flux[(k + 1) * 5 - 4] - flux[k * 5 - 4]) + 
                    dz2 * tz1 * (ue[(k + 1) * 5 - 4] - ue[
                    k * 5 - 4] * 2. + ue[(k - 1) * 5 - 4]);
	
	FRCT[3 + offset] = FRCT[3 + offset] + tz3 * .1 *
                    1. * (flux[(k + 1) * 5 - 3] - flux[k * 5 - 3]) + 
                    dz3 * tz1 * (ue[(k + 1) * 5 - 3] - ue[
                    k * 5 - 3] * 2. + ue[(k - 1) * 5 - 3]);
	
	FRCT[4 + offset] = FRCT[4 + offset] + tz3 * .1 *
                    1. * (flux[(k + 1) * 5 - 2] - flux[k * 5 - 2]) + 
                    dz4 * tz1 * (ue[(k + 1) * 5 - 2] - ue[
                    k * 5 - 2] * 2. + ue[(k - 1) * 5 - 2]);
	
	FRCT[5 + offset] = FRCT[5 + offset] + tz3 * .1 *
                    1. * (flux[(k + 1) * 5 - 1] - flux[k * 5 - 1]) + 
                    dz5 * tz1 * (ue[(k + 1) * 5 - 1] - ue[
                    k * 5 - 1] * 2. + ue[(k - 1) * 5 - 1]);
      }
      
/* ***fourth-order dissipation */
      
      for (m = 1; m <= 5; ++m) {
	
	FRCT[m + 5 * (x + y + Z_A[1]) - 1] -= dsspm * (ue[
                    m + 4] * 5. - ue[m + 9] * 4. + ue[m + 14]);
	
	FRCT[m + 5 * (x + y + Z_A[2]) - 1] -= dsspm * (ue[
                    m + 4] * -4. + ue[m + 9] * 6. - ue[m + 14] * 4. + ue[
                    m + 19]);
      }
      
      for (k = 4; k <= nz - 3; ++k) {
	ZO(z, k);
	offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m)
	  FRCT[m + offset] -= 
                    dsspm * (ue[m + (k - 2) * 5 - 6] - ue[m + (k - 1) 
                    * 5 - 6] * 4. + ue[m + k * 5 - 6] * 6. - ue[m + (
                    k + 1) * 5 - 6] * 4. + ue[m + (k + 2) * 5 - 6]);
      }
      for (m = 1; m <= 5; ++m) {
	
	FRCT[m + 5 * (x + y + Z_A[nz - 3]) - 1] -= dsspm * (ue[m + (nz - 4) * 5 - 6]
		    - ue[m + (nz - 3) * 5 - 6] * 4. + ue[m + (
                    nz - 2) * 5 - 6] * 6. - ue[m + (nz - 
                    1) * 5 - 6] * 4.);
	
	FRCT[m + 5 * (x + y + Z_A[nz - 2]) - 1] -= dsspm * (ue[m + (nz - 3) * 5 - 6] 
                    - ue[m + (nz - 2) * 5 - 6] * 4. + ue[m + (
                    nz - 1) * 5 - 6] * 5.);
      }
    }
  }
  return(0);
} /* erhs_ */

int rhs_(nx, ny, nz, tx, ty, tz, dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, dy5, dz1, dz2, dz3, dz4, dz5)
int nx, ny, nz;
double *tx, *ty, *tz, dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, dy5, dz1, dz2, dz3, dz4, dz5;
{
/* System generated locals */
  double d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8;
  
/* Local variables */
  double flux[5*S1], u21im1, u31im1, u41im1, u51im1, u21jm1;
  double u31jm1, u41jm1, u51jm1, u21km1, u31km1, u41km1, u51km1;
  int i, j, k, m;
  double q, u21, u31, u41, u21i, u31i, u41i, u51i, u21j, u31j;
  double u41j, u51j, u21k, u31k, u41k, u51k, tmp;
  double tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3;
  double dssp;

  int x, y, z;
  int offset;

/* ***compute the right hand sides */
  
  dssp = 1. / 4.;

  tx1 = tx[0];
  tx2 = tx[1];
  tx3 = tx[2];
  ty1 = ty[0];
  ty2 = ty[1];
  ty3 = ty[2];
  tz1 = tz[0];
  tz2 = tz[1];
  tz3 = tz[2];

  for (k = 1; k <= nz; ++k) {
    ZO(z, k);
    for (j = 1; j <= ny; ++j) {
      YO(y, j);
      for (i = 1; i <= nx; ++i) {
        XO(x, i);
	offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m)
	  RSD[m + offset] =  -FRCT[m + offset];
      }
    }
  }
  
/* ***xi-direction flux differences */
  
  for (k = 2; k <= nz - 1; ++k) {
    ZO(z, k)
    for (j = 2; j <= ny - 1; ++j) {
      YO(y, j)
      for (i = 1; i <= nx; ++i) {
        XO(x, i)
        offset = 5 * (x + y + z) - 1;
	flux[i * 5 - 5] = UA[2 + offset];
	u21 = UA[2 + offset] / UA[1 + offset];
	q = (UA[2 + offset] * UA[2 + offset] + UA[3 + offset] * UA[3 + offset] + UA[4 + offset] * UA[4 + offset]) * .5 / UA[1 + offset];
	
	flux[i * 5 - 4] = UA[2 + offset] 
                    * u21 + (UA[5 + offset] 
                    - q) * .4;
	flux[i * 5 - 3] = UA[3 + offset] * u21;
	flux[i * 5 - 2] = UA[4 + offset] * u21;
	flux[i * 5 - 1] = (UA[5 + offset]
                    * 1.4 - q * .4) * u21;
      }
      
      for (i = 2; i <= nx - 1; ++i) {
	XO(x, i)
	offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m)
	  RSD[m + offset] -= tx2 * (flux[m + (i + 1) * 5 - 6] - flux[m 
                    + (i - 1) * 5 - 6]);
      }
      
      for (i = 2; i <= nx; ++i) {
	XO(x, i)
        offset = 5 * (x + y + z) - 1;
	tmp = 1. / UA[1 + offset];
	
	u21i = tmp * UA[2 + offset];
	u31i = tmp * UA[3 + offset];
	u41i = tmp * UA[4 + offset];
	u51i = tmp * UA[5 + offset];
	
	tmp = 1. / UA[1 + 5 * (X_A[i-2] + y + z) - 1];
	u21im1 = tmp * UA[2 + 5 * (X_A[i-2] + y + z) - 1];
	u31im1 = tmp * UA[3 + 5 * (X_A[i-2] + y + z) - 1];
	u41im1 = tmp * UA[4 + 5 * (X_A[i-2] + y + z) - 1];
	u51im1 = tmp * UA[5 + 5 * (X_A[i-2] + y + z) - 1];
	
	flux[i * 5 - 4] = tx3 * 1.3333333333333333 * (u21i - u21im1);
	flux[i * 5 - 3] = tx3 * (u31i - u31im1);
	flux[i * 5 - 2] = tx3 * (u41i - u41im1);
	/* Computing 2nd power */
	d_1 = u21i;
	/* Computing 2nd power */
	d_2 = u31i;
	/* Computing 2nd power */
	d_3 = u41i;
	/* Computing 2nd power */
	d_4 = u21im1;
	/* Computing 2nd power */
	d_5 = u31im1;
	/* Computing 2nd power */
	d_6 = u41im1;
	/* Computing 2nd power */
	d_7 = u21i;
	/* Computing 2nd power */
	d_8 = u21im1;
	flux[i * 5 - 1] = tx3 * -.47999999999999987 * (d_1 * 
                    d_1 + d_2 * d_2 + d_3 * d_3 - (d_4 * d_4 + d_5 * d_5 
                    + d_6 * d_6)) + tx3 * .16666666666666666 * (
                    d_7 * d_7 - d_8 * d_8) + tx3 * 
                    1.9599999999999997 * (u51i - u51im1);
	
      }
      
      for (i = 2; i <= nx - 1; ++i) {
        XO(x, i)
        offset = 5 * (x + y + z) - 1;
	RSD[1 + offset] += dx1 * tx1 * (UA[1 + 5 * (X_A[i-2] + y + z) - 1]
		    - UA[1 + offset] * 2. + UA[1 + 5 * (X_A[i] + y + z) - 1]);
	
	RSD[2 + offset] = RSD[2 + offset] + tx3 * .1 * 
                    1. * (flux[(i + 1) * 5 - 4] - flux[i * 5 - 4]) + 
                    dx2 * tx1 * (UA[2 + 5 * (X_A[i-2] + y + z) - 1] - UA[2 + offset] * 2. + UA[2 + 5 * (X_A[i] + y + z) - 1]);
	
	RSD[3 + offset] = RSD[3 + offset] + tx3 * .1 * 
                    1. * (flux[(i + 1) * 5 - 3] - flux[i * 5 - 3]) + 
                    dx3 * tx1 * (UA[3 + 5 * (X_A[i-2] + y + z) - 1] - UA[3 + offset] * 2. + UA[3 + 5 * (X_A[i] + y + z) - 1]);
	
	RSD[4 + offset] = RSD[4 + offset] + tx3 * .1 * 
                    1. * (flux[(i + 1) * 5 - 2] - flux[i * 5 - 2]) + 
                    dx4 * tx1 * (UA[4 + 5 * (X_A[i-2] + y + z) - 1] - UA[4 + offset] * 2. + UA[4 + 5 * (X_A[i] + y + z) - 1]);
	
	RSD[5 + offset] = RSD[5 + offset] + tx3 * .1 * 
                    1. * (flux[(i + 1) * 5 - 1] - flux[i * 5 - 1]) + 
                    dx5 * tx1 * (UA[5 + 5 * (X_A[i-2] + y + z) - 1] - UA[5 + offset] * 2. + UA[5 + 5 * (X_A[i] + y + z) - 1]);
	
      }
      
/* ***Fourth-order dissipation */
      
      for (m = 1; m <= 5; ++m) {
	
	RSD[m + 5 * (X_A[1] + y + z) - 1] -= 
                    dssp * (UA[m + 5 * (X_A[1] + y + z) - 1] * 5. - UA[m + 5 * 
                    (X_A[2] + y + z) - 1] * 4. + UA[m + 5 * (X_A[3] + y + z) - 1]);
	
	RSD[m + 5 * (X_A[2] + y + z) - 1] -= dssp * (UA[m + 5 * (X_A[1] + y + z) - 1]
		  * -4. + UA[m + 5 * (X_A[2] + y + z) - 1] * 6. - UA[m + 5 * (X_A[3] +
                   y + z) - 1] * 4. + UA[m + 5 * (X_A[4] + y + z) - 1]);
      }
      
      for (i = 4; i <= nx - 3; ++i) {
        XO(x, i)
        offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m) {
	  
	  RSD[m + offset] -= dssp * (UA[m + 5 * (X_A[i-3] + y + z) - 1] - 
             UA[m + 5 * (X_A[i-2] + y + z) - 1] * 4. + UA[m + offset] * 6. - 
           UA[m + 5 * (X_A[i] + y + z) - 1] * 4. + UA[m + 5 * (X_A[i+1] + y + z) - 1]);
	}
      }
      
      for (m = 1; m <= 5; ++m) {
	
	RSD[m + 5 * (X_A[nx-3] + y + z) - 1] -= dssp * (UA[m + 5 * (X_A[nx-5] + y + 
            z) - 1] - UA[m + 5 * (X_A[nx-4] + y + z) - 1] * 4. + UA[m + 5 * 
            (X_A[nx-3] + y + z) - 1] * 6. - UA[m + 5 * (X_A[nx-2] + y + z) - 1] * 4.);
	
	RSD[m + 5 * (X_A[nx-2] + y + z) - 1] -= dssp * (UA[m + 5 * (X_A[nx-4] + y + z) - 1] - UA[m + 5 * (X_A[nx-3] + y + z) - 1] * 4. + UA[m + 5 * (X_A[nx-2] + y + z) - 1] * 5.);
      }
    }
  }
  
/* ***eta-direction flux differences */
  
  for (k = 2; k <= nz - 1; ++k) {
    ZO(z, k)
    for (i = 2; i <= nx - 1; ++i) {
      XO(x, i)
      for (j = 1; j <= ny; ++j) {
	YO(y, j)
        offset = 5 * (x + y + z) - 1;
	flux[j * 5 - 5] = UA[3 + offset];
	u31 = UA[3 + offset] / UA[1 + offset];
                    
	q = (UA[2 + offset] * UA[2 + offset] + UA[3 + offset] * UA[3 + offset] + UA[4 + offset] * UA[4 + offset]) * .5 / UA[1 + offset];
	
	flux[j * 5 - 4] = UA[2 + offset] * u31;
	
	flux[j * 5 - 3] = UA[3 + offset] 
                    * u31 + (UA[5 + offset] 
		   - q) * .4;
                    
	flux[j * 5 - 2] = UA[4 + offset] * u31;
	
	flux[j * 5 - 1] = (UA[5 + offset]
                    * 1.4 - q * .4) * u31;
      }
      
      for (j = 2; j <= ny - 1; ++j) {
	YO(y, j)
        offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m)
	  RSD[m + offset] -= ty2 * (flux[m + (j + 1) * 5 - 6] - 
                             flux[m + (j - 1) * 5 - 6]);
      }
      
      for (j = 2; j <= ny; ++j) {
	YO(y, j)
        offset = 5 * (x + y + z) - 1;
	tmp = 1. / UA[1 + offset];
	
	u21j = tmp * UA[2 + offset];
	u31j = tmp * UA[3 + offset];
	u41j = tmp * UA[4 + offset];
	u51j = tmp * UA[5 + offset];
	
	tmp = 1. / UA[1 + 5 * (x + Y_A[j-2] + z) - 1];
	
	u21jm1 = tmp * UA[2 + 5 * (x + Y_A[j-2] + z) - 1];
	u31jm1 = tmp * UA[3 + 5 * (x + Y_A[j-2] + z) - 1];
	u41jm1 = tmp * UA[4 + 5 * (x + Y_A[j-2] + z) - 1];
	u51jm1 = tmp * UA[5 + 5 * (x + Y_A[j-2] + z) - 1];
	
	flux[j * 5 - 4] = ty3 * (u21j - u21jm1);
	flux[j * 5 - 3] = ty3 * 1.3333333333333333 * (u31j - u31jm1);
	flux[j * 5 - 2] = ty3 * (u41j - u41jm1);
	/* Computing 2nd power */
	d_1 = u21j;
	/* Computing 2nd power */
	d_2 = u31j;
	/* Computing 2nd power */
	d_3 = u41j;
	/* Computing 2nd power */
	d_4 = u21jm1;
	/* Computing 2nd power */
	d_5 = u31jm1;
	/* Computing 2nd power */
	d_6 = u41jm1;
	/* Computing 2nd power */
	d_7 = u31j;
	/* Computing 2nd power */
	d_8 = u31jm1;
	flux[j * 5 - 1] = ty3 * -.47999999999999987 * (d_1 * 
                    d_1 + d_2 * d_2 + d_3 * d_3 - (d_4 * d_4 + d_5 * d_5 
                    + d_6 * d_6)) + ty3 * .16666666666666666 * (
                    d_7 * d_7 - d_8 * d_8) + ty3 * 
                    1.9599999999999997 * (u51j - u51jm1);
	
      }
      
      for (j = 2; j <= ny - 1; ++j) {
	YO(y, j)
        offset = 5 * (x + y + z) - 1;
	RSD[1 + offset] += dy1 * 
                    ty1 * (UA[1 + 5 * (x + Y_A[j-2] + z) - 1] - UA[1 + offset]
                    * 2. + UA[1 + 5 * (x + Y_A[j] + z) - 1]);
	
	RSD[2 + offset] = RSD[2 + offset] + ty3 * .1 * 
                    1. * (flux[(j + 1) * 5 - 4] - flux[j * 5 - 4]) + 
                    dy2 * ty1 * (UA[2 + 5 * (x + Y_A[j-2] + z) - 1] - UA[2 + offset] * 2. + UA[2 + 5 * (x + Y_A[j] + z) - 1]);
	
	RSD[3 + offset] = RSD[3 + offset] + ty3 * .1 * 
                    1. * (flux[(j + 1) * 5 - 3] - flux[j * 5 - 3]) + 
                    dy3 * ty1 * (UA[3 + 5 * (x + Y_A[j-2] + z) - 1] - UA[3 + offset] * 2. + UA[3 + 5 * (x + Y_A[j] + z) - 1]);
	
	RSD[4 + offset] = RSD[4 + offset] + ty3 * .1 * 
                    1. * (flux[(j + 1) * 5 - 2] - flux[j * 5 - 2]) + 
                    dy4 * ty1 * (UA[4 + 5 * (x + Y_A[j-2] + z) - 1] - UA[4 + offset] * 2. + UA[4 + 5 * (x + Y_A[j] + z) - 1]);
	
	RSD[5 + offset] = RSD[5 + offset] + ty3 * .1 * 
                    1. * (flux[(j + 1) * 5 - 1] - flux[j * 5 - 1]) + 
                    dy5 * ty1 * (UA[5 + 5 * (x + Y_A[j-2] + z) - 1] - UA[5 + offset] * 2. + UA[5 + 5 * (x + Y_A[j] + z) - 1]);
      }
      
/* ***fourth-order dissipation */
      
      for (m = 1; m <= 5; ++m) {
	
	RSD[m + 5 * (x + Y_A[1] + z) - 1] -= 
                    dssp * (UA[m + 5 * (x + Y_A[1] + z) - 1] * 5. - UA[m + 5 * (x + Y_A[2] + z) - 1] * 4. + UA[m + 5 * (x + Y_A[3] + z) - 1]);
	
	RSD[m + 5 * (x + Y_A[2] + z) - 1] -= 
                    dssp * (UA[m + 5 * (x + Y_A[1] + z) - 1] * -4. + UA[m + 5 * (x + Y_A[2] + z) - 1] * 6. - UA[m + 5 * (x + Y_A[3] + z) - 1] * 4. + UA[m + 5 * (x + Y_A[4] + z) - 1]);
      }
      
      for (j = 4; j <= ny - 3; ++j) {
        YO(y, j)
        offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m)
	  RSD[m + offset] -= 
                    dssp * (UA[m + 5 * (x + Y_A[j-3] + z) - 1] - UA[m + 5 * (x + Y_A[j-2] + z) - 1] * 4. + UA[m + offset] * 6. - UA[m + 5 * (x + Y_A[j] + z) - 1] * 4. + 
                    UA[m + 5 * (x + Y_A[j+1] + z) - 1]);
      }
      
      for (m = 1; m <= 5; ++m) {
	RSD[m + 5 * (x + Y_A[ny-3] + z) - 1]
                    -= dssp * (UA[m + 5 * (x + Y_A[ny-5] + z) - 1] - UA[m + 5 * (x + Y_A[ny-4] + z) - 1] * 4. + 
                    UA[m + 5 * (x + Y_A[ny-3] + z) - 1] * 6. - UA[m + 5 * (x + Y_A[ny-2] + z) - 1] * 4.);
	
	RSD[m + 5 * (x + Y_A[ny-2] + z) - 1]
                    -= dssp * (UA[m + 5 * (x + Y_A[ny-4] + z) - 1] - UA[m + 5 * (x + Y_A[ny-3] + z) - 1] * 4. + 
                    UA[m + 5 * (x + Y_A[ny-2] + z) - 1] * 5.);
      }
    }
  }
  
/* ***zeta-direction flux differences */
  
  for (j = 2; j <= ny - 1; ++j) {
    YO(y, j)
    for (i = 2; i <= nx - 1; ++i) {
      XO(x, i)
      for (k = 1; k <= nz; ++k) {
        ZO(z, k)
        offset = 5 * (x + y + z) - 1;
	
	flux[k * 5 - 5] = UA[4 + offset];
	u41 = UA[4 + offset] / UA[1 + offset];
	
	q = (UA[2 + offset] * UA[2 + offset] + UA[3 + offset] * UA[3 + offset] + UA[4 + offset] * UA[4 + offset]) * .5 / UA[1 + offset];
	
	flux[k * 5 - 4] = UA[2 + offset] * u41;
	flux[k * 5 - 3] = UA[3 + offset] * u41;
	flux[k * 5 - 2] = UA[4 + offset] 
                    * u41 + (UA[5 + offset] 
                    - q) * .4;
	
	flux[k * 5 - 1] = (UA[5 + offset]
                    * 1.4 - q * .4) * u41;
      }
      
      for (k = 2; k <= nz - 1; ++k) {
	ZO(z, k)
        offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m)
	  RSD[m + offset] -= tz2 * (flux[m + (k + 1) * 5 - 6] - flux[m 
                    + (k - 1) * 5 - 6]);
      }
      
      for (k = 2; k <= nz; ++k) {
	ZO(z, k)
        offset = 5 * (x + y + z) - 1;
	tmp = 1. / UA[1 + offset];
	
	u21k = tmp * UA[2 + offset];
	u31k = tmp * UA[3 + offset];
	u41k = tmp * UA[4 + offset];
	u51k = tmp * UA[5 + offset];
	
	tmp = 1. / UA[1 + 5 * (x + y + Z_A[k-2]) - 1];
	
	u21km1 = tmp * UA[2 + 5 * (x + y + Z_A[k-2]) - 1];
	u31km1 = tmp * UA[3 + 5 * (x + y + Z_A[k-2]) - 1];
	u41km1 = tmp * UA[4 + 5 * (x + y + Z_A[k-2]) - 1];
	u51km1 = tmp * UA[5 + 5 * (x + y + Z_A[k-2]) - 1];
	
	flux[k * 5 - 4] = tz3 * (u21k - u21km1);
	flux[k * 5 - 3] = tz3 * (u31k - u31km1);
	flux[k * 5 - 2] = tz3 * 1.3333333333333333 * (u41k - u41km1);
	/* Computing 2nd power */
	d_1 = u21k;
	/* Computing 2nd power */
	d_2 = u31k;
	/* Computing 2nd power */
	d_3 = u41k;
	/* Computing 2nd power */
	d_4 = u21km1;
	/* Computing 2nd power */
	d_5 = u31km1;
	/* Computing 2nd power */
	d_6 = u41km1;
	/* Computing 2nd power */
	d_7 = u41k;
	/* Computing 2nd power */
	d_8 = u41km1;
	flux[k * 5 - 1] = tz3 * -.47999999999999987 * (d_1 * 
                    d_1 + d_2 * d_2 + d_3 * d_3 - (d_4 * d_4 + d_5 * d_5 
                    + d_6 * d_6)) + tz3 * .16666666666666666 * (
                    d_7 * d_7 - d_8 * d_8) + tz3 * 
                    1.9599999999999997 * (u51k - u51km1);
      }
      
      for (k = 2; k <= nz - 1; ++k) {
	ZO(z, k)
        offset = 5 * (x + y + z) - 1;
	RSD[1 + offset] += dz1 * 
                    tz1 * (UA[1 + 5 * (x + y + Z_A[k-2]) - 1] - UA[1 + offset] * 2. + UA[1 + 5 * (x + y + Z_A[k]) - 1]);
	
	RSD[2 + offset] = RSD[2 + offset] + tz3 * .1 * 
                    1. * (flux[(k + 1) * 5 - 4] - flux[k * 5 - 4]) + 
                    dz2 * tz1 * (UA[2 + 5 * (x + y + Z_A[k-2]) - 1] - UA[2 + offset] * 2. + UA[2 + 5 * (x + y + Z_A[k]) - 1]);
	
	RSD[3 + offset] = RSD[3 + offset] + tz3 * .1 * 
                    1. * (flux[(k + 1) * 5 - 3] - flux[k * 5 - 3]) + 
                    dz3 * tz1 * (UA[3 + 5 * (x + y + Z_A[k-2]) - 1] - UA[3 + offset] * 2. + UA[3 + 5 * (x + y + Z_A[k]) - 1]);
	
	RSD[4 + offset] = RSD[4 + offset] + tz3 * .1 * 
                    1. * (flux[(k + 1) * 5 - 2] - flux[k * 5 - 2]) + 
                    dz4 * tz1 * (UA[4 + 5 * (x + y + Z_A[k-2]) - 1] - UA[4 + offset] * 2. + UA[4 + 5 * (x + y + Z_A[k]) - 1]);
	
	RSD[5 + offset] = RSD[5 + offset] + tz3 * .1 * 
                    1. * (flux[(k + 1) * 5 - 1] - flux[k * 5 - 1]) + 
                    dz5 * tz1 * (UA[5 + 5 * (x + y + Z_A[k-2]) - 1] - UA[5 + offset] * 2. + UA[5 + 5 * (x + y + Z_A[k]) - 1]);
      }
      
/* ***fourth-order dissipation */
      
      for (m = 1; m <= 5; ++m) {
	
	RSD[m + 5 * (x + y + Z_A[1]) - 1] -= dssp *
                    (UA[m + 5 * (x + y + Z_A[1]) - 1] * 5. - 
                    UA[m + 5 * (x + y + Z_A[2]) - 1] * 4. + 
                    UA[m + 5 * (x + y + Z_A[3]) - 1]);
	
	RSD[m + 5 * (x + y + Z_A[2]) - 1] -= dssp *
                    (UA[m + 5 * (x + y + Z_A[1]) - 1] * -4. + 
                    UA[m + 5 * (x + y + Z_A[2]) - 1] * 6. - 
                    UA[m + 5 * (x + y + Z_A[3]) - 1] * 4. + 
                    UA[m + 5 * (x + y + Z_A[4]) - 1]);
      }
      
      for (k = 4; k <= nz - 3; ++k) {
        ZO(z, k)
        offset = 5 * (x + y + z) - 1;
	for (m = 1; m <= 5; ++m)
	  RSD[m + offset] -= dssp * (UA[m + 5 * (x + y + Z_A[k-3]) - 1] - 
                    UA[m + 5 * (x + y + Z_A[k-2]) - 1] * 4. + UA[m + offset] * 6. - 
                    UA[m + 5 * (x + y + Z_A[k]) - 1] * 4. 
                    + UA[m + 5 * (x + y + Z_A[k+1]) - 1]);
      }
      
      for (m = 1; m <= 5; ++m) {
	RSD[m + 5 * (x + y + Z_A[nz-3]) - 1] -= dssp * (UA[m + 5 * (x + y + 
             Z_A[nz-5]) - 1] - UA[m + 5 * (x + y + Z_A[nz-4]) - 1] * 4. + UA[m + 5 * 
            (x + y + Z_A[nz-3]) - 1] * 6. - UA[m + 5 * (x + y + Z_A[nz-2]) - 1] * 4.);
	
	RSD[m + 5 * (x + y + Z_A[nz-2]) - 1] -= dssp * (UA[m + 5 * (x + y + 
                    Z_A[nz-4]) - 1] - UA[m + 5 * (x + y + Z_A[nz-3]) - 1] * 
                    4. + UA[m + 5 * (x + y + Z_A[nz-2]) - 1] * 5.);
      }
    }
  }
  
  return(0);
} /* rhs_ */

