/* --- jac.c ---
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

int jacx_(dt, tx1, tx2, dx1, dx2, dx3, dx4, dx5, nx, ny, nz)
double dt, tx1, tx2, dx1, dx2, dx3, dx4, dx5;
int nx, ny, nz;
{
/* System generated locals */
  double d_1, d_2, d_3;
  
/* Local variables */
  double fjac[5*5*S1];
  double njac[5*5*S1];
  int i, j, k, l, m;
  double c34, r43, c1345, tmp1, tmp2, tmp3;
  int x, y, z;
  int offset;
  
/* ***form the block tridiagonal system for xi-direction sweep */
  
  r43 = 1.3333333333333333;
  c34 = .10000000000000001;
  c1345 = .19599999999999998;
  
  for (k = 2; k <= nz - 1; ++k) {
    ZO(z, k)
    for (j = 2; j <= ny - 1; ++j) {
      YO(y, j)      
/* ***compute the xi-direction flux jacobian */
      
      for (i = 1; i <= nx; ++i) {
	XO(x, i)
        offset = 5 * (x + y + z) - 1;
	fjac[(i * 5 + 1) * 5 - 30] = 0.;
	fjac[(i * 5 + 2) * 5 - 30] = 1.;
	fjac[(i * 5 + 3) * 5 - 30] = 0.;
	fjac[(i * 5 + 4) * 5 - 30] = 0.;
	fjac[(i * 5 + 5) * 5 - 30] = 0.;
	
	/* Computing 2nd power */
	d_1 = UA[2 + offset] / UA[1 + offset];
	fjac[(i * 5 + 1) * 5 - 29] = -(d_1 * d_1) + (UA[2 + offset] * UA[2 + offset] + UA[3 + offset] * UA[3 + offset] + UA[4 + offset] *
                    UA[4 + offset]) * 
                    .20000000000000001 / (UA[1 + offset] * UA[1 + offset]);
	fjac[(i * 5 + 2) * 5 - 29] = UA[2 + offset] / UA[1 + offset] * 1.6000000000000001;
	fjac[(i * 5 + 3) * 5 - 29] = UA[3 + offset] / UA[1 + offset] * -.4;
	fjac[(i * 5 + 4) * 5 - 29] = UA[4 + offset] / UA[1 + offset] * -.4;
	fjac[(i * 5 + 5) * 5 - 29] = .4;
	
	fjac[(i * 5 + 1) * 5 - 28] = -(UA[2 + offset] * UA[3 + offset]) / (UA[1 + offset]
                    * UA[1 + offset]);
	fjac[(i * 5 + 2) * 5 - 28] = UA[3 + offset] / UA[1 + offset];
	fjac[(i * 5 + 3) * 5 - 28] = UA[2 + offset] / UA[1 + offset];
	fjac[(i * 5 + 4) * 5 - 28] = 0.;
	fjac[(i * 5 + 5) * 5 - 28] = 0.;
	
	fjac[(i * 5 + 1) * 5 - 27] = -(UA[2 + offset] * UA[4 + offset]) / (UA[1 + offset]
                    * UA[1 + offset]);
	fjac[(i * 5 + 2) * 5 - 27] = UA[4 + offset] / UA[1 + offset];
	fjac[(i * 5 + 3) * 5 - 27] = 0.;
	fjac[(i * 5 + 4) * 5 - 27] = UA[2 + offset] / UA[1 + offset];
	fjac[(i * 5 + 5) * 5 - 27] = 0.;
	
	fjac[(i * 5 + 1) * 5 - 26] = ((UA[2 + offset] * UA[2 + offset] + UA[3 + offset] *
                    UA[3 + offset] + 
                    UA[4 + offset] * 
                    UA[4 + offset]) * .4 / (
                    UA[1 + offset] * 
                    UA[1 + offset]) - 
                    UA[5 + offset] / 
                    UA[1 + offset] * 1.4) * (
                    UA[2 + offset] / 
                    UA[1 + offset]);
	fjac[(i * 5 + 2) * 5 - 26] = UA[5 + offset] / UA[1 + offset] * 1.4 - (UA[2 + offset] * 3. * UA[2 + offset] + UA[3 + offset] * 
                    UA[3 + offset] + 
                    UA[4 + offset] * 
                    UA[4 + offset]) / (
                    UA[1 + offset] * 
                    UA[1 + offset]) * 
                    .20000000000000001;
	fjac[(i * 5 + 3) * 5 - 26] = UA[3 + offset] * UA[2 + offset] * -.4 / (UA[1 + offset] * UA[1 + offset]);
	fjac[(i * 5 + 4) * 5 - 26] = UA[4 + offset] * UA[2 + offset] * -.4 / (UA[1 + offset] * UA[1 + offset]);
	fjac[(i * 5 + 5) * 5 - 26] = UA[2 + offset] / UA[1 + offset] * 1.4;
	
	tmp1 = 1. / UA[1 + offset];
	tmp2 = tmp1 * tmp1;
	tmp3 = tmp1 * tmp2;
	
	njac[(i * 5 + 1) * 5 - 30] = 0.;
	njac[(i * 5 + 2) * 5 - 30] = 0.;
	njac[(i * 5 + 3) * 5 - 30] = 0.;
	njac[(i * 5 + 4) * 5 - 30] = 0.;
	njac[(i * 5 + 5) * 5 - 30] = 0.;
	
	njac[(i * 5 + 1) * 5 - 29] = -r43 * c34 * tmp2 * UA[2 + offset];
	njac[(i * 5 + 2) * 5 - 29] = r43 * c34 * tmp1;
	njac[(i * 5 + 3) * 5 - 29] = 0.;
	njac[(i * 5 + 4) * 5 - 29] = 0.;
	njac[(i * 5 + 5) * 5 - 29] = 0.;
	
	njac[(i * 5 + 1) * 5 - 28] = -c34 * tmp2 * UA[3 + offset];
	njac[(i * 5 + 2) * 5 - 28] = 0.;
	njac[(i * 5 + 3) * 5 - 28] = c34 * tmp1;
	njac[(i * 5 + 4) * 5 - 28] = 0.;
	njac[(i * 5 + 5) * 5 - 28] = 0.;
	
	njac[(i * 5 + 1) * 5 - 27] = -c34 * tmp2 * UA[4 + offset];
	njac[(i * 5 + 2) * 5 - 27] = 0.;
	njac[(i * 5 + 3) * 5 - 27] = 0.;
	njac[(i * 5 + 4) * 5 - 27] = c34 * tmp1;
	njac[(i * 5 + 5) * 5 - 27] = 0.;
	
	/* Computing 2nd power */
	d_1 = UA[2 + offset];
	/* Computing 2nd power */
	d_2 = UA[3 + offset];
	/* Computing 2nd power */
	d_3 = UA[4 + offset];
	njac[(i * 5 + 1) * 5 - 26] = -(r43 * c34 - c1345) * tmp3 * (
                    d_1 * d_1) - (c34 - c1345) * tmp3 * (d_2 * d_2) - (
                    c34 - c1345) * tmp3 * (d_3 * d_3) - c1345 * tmp2 * 
                    UA[5 + offset];
	
	njac[(i * 5 + 2) * 5 - 26] = (r43 * c34 - c1345) * tmp2 * 
                    UA[2 + offset];
	njac[(i * 5 + 3) * 5 - 26] = (c34 - c1345) * tmp2 * UA[3 + offset];
	njac[(i * 5 + 4) * 5 - 26] = (c34 - c1345) * tmp2 * UA[4 + offset];
	njac[(i * 5 + 5) * 5 - 26] = c1345 * tmp1;
      }
      
/* ***dirichlet boundary conditions */
      
      for (l = 1; l <= 5; ++l) {
	for (m = 1; m <= 5; ++m) {
	  AA(m,l,(5 * (X_A[0]+y+z))) = 0.;
	  BA(m,l,(5 * (X_A[0]+y+z))) = 0.;
	  CA(m,l,(5 * (X_A[0]+y+z))) = 0.;
	}
	BA(l, l, (5 * (X_A[0]+y+z))) = 1.;
      }

      for (i = 2; i <= nx - 1; ++i) {
	XO(x, i)
	offset = 5 * (x + y + z);
	tmp1 = dt * tx1;
	tmp2 = dt * tx2;
	
	AA(1,1,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 1) * 5 - 30] - tmp1 * 
                    njac[((i - 1) * 5 + 1) * 5 - 30] - tmp1 * dx1;
	AA(1,2,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 2) * 5 - 30] - tmp1 * 
                    njac[((i - 1) * 5 + 2) * 5 - 30];
	AA(1,3,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 3) * 5 - 30] - tmp1 * 
                    njac[((i - 1) * 5 + 3) * 5 - 30];
	AA(1,4,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 4) * 5 - 30] - tmp1 * 
                    njac[((i - 1) * 5 + 4) * 5 - 30];
	AA(1,5,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 5) * 5 - 30] - tmp1 * 
                    njac[((i - 1) * 5 + 5) * 5 - 30];
	
	AA(2,1,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 1) * 5 - 29] - tmp1 * 
                    njac[((i - 1) * 5 + 1) * 5 - 29];
	AA(2,2,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 2) * 5 - 29] - tmp1 * 
                    njac[((i - 1) * 5 + 2) * 5 - 29] - tmp1 * dx2;
	AA(2,3,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 3) * 5 - 29] - tmp1 * 
                    njac[((i - 1) * 5 + 3) * 5 - 29];
	AA(2,4,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 4) * 5 - 29] - tmp1 * 
                    njac[((i - 1) * 5 + 4) * 5 - 29];
	AA(2,5,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 5) * 5 - 29] - tmp1 * 
                    njac[((i - 1) * 5 + 5) * 5 - 29];
	
	AA(3,1,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 1) * 5 - 28] - tmp1 * 
                    njac[((i - 1) * 5 + 1) * 5 - 28];
	AA(3,2,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 2) * 5 - 28] - tmp1 * 
                    njac[((i - 1) * 5 + 2) * 5 - 28];
	AA(3,3,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 3) * 5 - 28] - tmp1 * 
                    njac[((i - 1) * 5 + 3) * 5 - 28] - tmp1 * dx3;
	AA(3,4,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 4) * 5 - 28] - tmp1 * 
                    njac[((i - 1) * 5 + 4) * 5 - 28];
	AA(3,5,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 5) * 5 - 28] - tmp1 * 
                    njac[((i - 1) * 5 + 5) * 5 - 28];
	
	AA(4,1,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 1) * 5 - 27] - tmp1 * 
                    njac[((i - 1) * 5 + 1) * 5 - 27];
	AA(4,2,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 2) * 5 - 27] - tmp1 * 
                    njac[((i - 1) * 5 + 2) * 5 - 27];
	AA(4,3,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 3) * 5 - 27] - tmp1 * 
                    njac[((i - 1) * 5 + 3) * 5 - 27];
	AA(4,4,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 4) * 5 - 27] - tmp1 * 
                    njac[((i - 1) * 5 + 4) * 5 - 27] - tmp1 * dx4;
	AA(4,5,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 5) * 5 - 27] - tmp1 * 
                    njac[((i - 1) * 5 + 5) * 5 - 27];
	
	AA(5,1,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 1) * 5 - 26] - tmp1 * 
                    njac[((i - 1) * 5 + 1) * 5 - 26];
	AA(5,2,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 2) * 5 - 26] - tmp1 * 
                    njac[((i - 1) * 5 + 2) * 5 - 26];
	AA(5,3,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 3) * 5 - 26] - tmp1 * 
                    njac[((i - 1) * 5 + 3) * 5 - 26];
	AA(5,4,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 4) * 5 - 26] - tmp1 * 
                    njac[((i - 1) * 5 + 4) * 5 - 26];
	AA(5,5,offset) = 
                    -tmp2 * fjac[((i - 1) * 5 + 5) * 5 - 26] - tmp1 * 
                    njac[((i - 1) * 5 + 5) * 5 - 26] - tmp1 * dx5;
	
	BA(1,1,offset) = tmp1 
                    * 2. * njac[(i * 5 + 1) * 5 - 30] + 1. + tmp1 * 2. * 
                    dx1;
	BA(1,2,offset) = tmp1 
                    * 2. * njac[(i * 5 + 2) * 5 - 30];
	BA(1,3,offset) = tmp1 
                    * 2. * njac[(i * 5 + 3) * 5 - 30];
	BA(1,4,offset) = tmp1 
                    * 2. * njac[(i * 5 + 4) * 5 - 30];
	BA(1,5,offset) = tmp1 
                    * 2. * njac[(i * 5 + 5) * 5 - 30];
	
	BA(2,1,offset) = tmp1 
                    * 2. * njac[(i * 5 + 1) * 5 - 29];
	BA(2,2,offset) = tmp1 
                    * 2. * njac[(i * 5 + 2) * 5 - 29] + 1. + tmp1 * 2. * 
                    dx2;
	BA(2,3,offset) = tmp1 
                    * 2. * njac[(i * 5 + 3) * 5 - 29];
	BA(2,4,offset) = tmp1 
                    * 2. * njac[(i * 5 + 4) * 5 - 29];
	BA(2,5,offset) = tmp1 
                    * 2. * njac[(i * 5 + 5) * 5 - 29];
	
	BA(3,1,offset) = tmp1 
                    * 2. * njac[(i * 5 + 1) * 5 - 28];
	BA(3,2,offset) = tmp1 
                    * 2. * njac[(i * 5 + 2) * 5 - 28];
	BA(3,3,offset) = tmp1 
                    * 2. * njac[(i * 5 + 3) * 5 - 28] + 1. + tmp1 * 2. * 
                    dx3;
	BA(3,4,offset) = tmp1 
                    * 2. * njac[(i * 5 + 4) * 5 - 28];
	BA(3,5,offset) = tmp1 
                    * 2. * njac[(i * 5 + 5) * 5 - 28];
	
	BA(4,1,offset) = tmp1 
                    * 2. * njac[(i * 5 + 1) * 5 - 27];
	BA(4,2,offset) = tmp1 
                    * 2. * njac[(i * 5 + 2) * 5 - 27];
	BA(4,3,offset) = tmp1 
                    * 2. * njac[(i * 5 + 3) * 5 - 27];
	BA(4,4,offset) = tmp1 
                    * 2. * njac[(i * 5 + 4) * 5 - 27] + 1. + tmp1 * 2. * 
                    dx4;
	BA(4,5,offset) = tmp1 
                    * 2. * njac[(i * 5 + 5) * 5 - 27];
	
	BA(5,1,offset) = tmp1 
                    * 2. * njac[(i * 5 + 1) * 5 - 26];
	BA(5,2,offset) = tmp1 
                    * 2. * njac[(i * 5 + 2) * 5 - 26];
	BA(5,3,offset) = tmp1 
                    * 2. * njac[(i * 5 + 3) * 5 - 26];
	BA(5,4,offset) = tmp1 
                    * 2. * njac[(i * 5 + 4) * 5 - 26];
	BA(5,5,offset) = tmp1 
                    * 2. * njac[(i * 5 + 5) * 5 - 26] + 1. + tmp1 * 2. * 
                    dx5;
	
	CA(1,1,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 1) * 5 - 30] - tmp1 * njac[((i 
                    + 1) * 5 + 1) * 5 - 30] - tmp1 * dx1;
	CA(1,2,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 2) * 5 - 30] - tmp1 * njac[((i 
                    + 1) * 5 + 2) * 5 - 30];
	CA(1,3,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 3) * 5 - 30] - tmp1 * njac[((i 
                    + 1) * 5 + 3) * 5 - 30];
	CA(1,4,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 4) * 5 - 30] - tmp1 * njac[((i 
                    + 1) * 5 + 4) * 5 - 30];
	CA(1,5,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 5) * 5 - 30] - tmp1 * njac[((i 
                    + 1) * 5 + 5) * 5 - 30];
	
	CA(2,1,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 1) * 5 - 29] - tmp1 * njac[((i 
                    + 1) * 5 + 1) * 5 - 29];
	CA(2,2,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 2) * 5 - 29] - tmp1 * njac[((i 
                    + 1) * 5 + 2) * 5 - 29] - tmp1 * dx2;
	CA(2,3,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 3) * 5 - 29] - tmp1 * njac[((i 
                    + 1) * 5 + 3) * 5 - 29];
	CA(2,4,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 4) * 5 - 29] - tmp1 * njac[((i 
                    + 1) * 5 + 4) * 5 - 29];
	CA(2,5,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 5) * 5 - 29] - tmp1 * njac[((i 
                    + 1) * 5 + 5) * 5 - 29];
	
	CA(3,1,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 1) * 5 - 28] - tmp1 * njac[((i 
                    + 1) * 5 + 1) * 5 - 28];
	CA(3,2,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 2) * 5 - 28] - tmp1 * njac[((i 
                    + 1) * 5 + 2) * 5 - 28];
	CA(3,3,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 3) * 5 - 28] - tmp1 * njac[((i 
                    + 1) * 5 + 3) * 5 - 28] - tmp1 * dx3;
	CA(3,4,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 4) * 5 - 28] - tmp1 * njac[((i 
                    + 1) * 5 + 4) * 5 - 28];
	CA(3,5,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 5) * 5 - 28] - tmp1 * njac[((i 
                    + 1) * 5 + 5) * 5 - 28];
	
	CA(4,1,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 1) * 5 - 27] - tmp1 * njac[((i 
                    + 1) * 5 + 1) * 5 - 27];
	CA(4,2,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 2) * 5 - 27] - tmp1 * njac[((i 
                    + 1) * 5 + 2) * 5 - 27];
	CA(4,3,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 3) * 5 - 27] - tmp1 * njac[((i 
                    + 1) * 5 + 3) * 5 - 27];
	CA(4,4,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 4) * 5 - 27] - tmp1 * njac[((i 
                    + 1) * 5 + 4) * 5 - 27] - tmp1 * dx4;
	CA(4,5,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 5) * 5 - 27] - tmp1 * njac[((i 
                    + 1) * 5 + 5) * 5 - 27];
	
	CA(5,1,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 1) * 5 - 26] - tmp1 * njac[((i 
                    + 1) * 5 + 1) * 5 - 26];
	CA(5,2,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 2) * 5 - 26] - tmp1 * njac[((i 
                    + 1) * 5 + 2) * 5 - 26];
	CA(5,3,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 3) * 5 - 26] - tmp1 * njac[((i 
                    + 1) * 5 + 3) * 5 - 26];
	CA(5,4,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 4) * 5 - 26] - tmp1 * njac[((i 
                    + 1) * 5 + 4) * 5 - 26];
	CA(5,5,offset) = tmp2 
                    * fjac[((i + 1) * 5 + 5) * 5 - 26] - tmp1 * njac[((i 
                    + 1) * 5 + 5) * 5 - 26] - tmp1 * dx5;
      }
      
/* ***dirichlet boundary conditions */
      
      for (l = 1; l <= 5; ++l) {
	for (m = 1; m <= 5; ++m) {
	  AA(m,l,(5*(X_A[nx-1]+y+z))) = 0.;
	  BA(m,l,(5*(X_A[nx-1]+y+z))) = 0.;
	  CA(m,l,(5*(X_A[nx-1]+y+z))) = 0.;
	}
	BA(l,l,(5*(X_A[nx-1]+y+z))) = 1.;
      }
    }
  }
  
  return(0);
} /* jacx_ */


int jacy_(dt, ty1, ty2, dy1, dy2, dy3, dy4, dy5, nx, ny, nz)
double dt, ty1, ty2, dy1, dy2, dy3, dy4, dy5;
int nx, ny, nz;
{
/* System generated locals */
  double d_1, d_2, d_3;
  
/* Local variables */
  double fjac[5*5*S2];
  double qjac[5*5*S2];
  int i, j, k, l, m;
  double c34, r43, c1345, tmp1, tmp2, tmp3;
  int x, y, z, offset;
  
/* ***form the block tridiagonal system for eta-direction sweep */
  
  r43 = 1.3333333333333333;
  c34 = .10000000000000001;
  c1345 = .19599999999999998;
  
  for (k = 2; k <= nz - 1; ++k) {
    ZO(z, k)
    for (i = 2; i <= nx - 1; ++i) {
      XO(x, i)
/* ***compute the eta-direction flux jacobians */
      
      for (j = 1; j <= ny; ++j) {
	YO(y, j)
        offset = 5 * (x + y + z) - 1;
	fjac[(j * 5 + 1) * 5 - 30] = 0.;
	fjac[(j * 5 + 2) * 5 - 30] = 0.;
	fjac[(j * 5 + 3) * 5 - 30] = 1.;
	fjac[(j * 5 + 4) * 5 - 30] = 0.;
	fjac[(j * 5 + 5) * 5 - 30] = 0.;
	
	fjac[(j * 5 + 1) * 5 - 29] = -(UA[2 + offset] * UA[3 + offset]) / (UA[1 + offset]
                    * UA[1 + offset]);
	fjac[(j * 5 + 2) * 5 - 29] = UA[3 + offset] / UA[1 + offset];
	fjac[(j * 5 + 3) * 5 - 29] = UA[2 + offset] / UA[1 + offset];
	fjac[(j * 5 + 4) * 5 - 29] = 0.;
	fjac[(j * 5 + 5) * 5 - 29] = 0.;
	
	/* Computing 2nd power */
	d_1 = UA[3 + offset] / UA[1 + offset];
	fjac[(j * 5 + 1) * 5 - 28] = -(d_1 * d_1) + (UA[2 + offset] * UA[2 + offset] + UA[3 + offset] * UA[3 + offset] + UA[4 + offset] *
                    UA[4 + offset]) / (
                    UA[1 + offset] * 
                    UA[1 + offset]) * 
                    .20000000000000001;
	fjac[(j * 5 + 2) * 5 - 28] = UA[2 + offset] / UA[1 + offset] * -.4;
	fjac[(j * 5 + 3) * 5 - 28] = UA[3 + offset] / UA[1 + offset] * 1.6000000000000001;
	fjac[(j * 5 + 4) * 5 - 28] = UA[4 + offset] / UA[1 + offset] * -.4;
	fjac[(j * 5 + 5) * 5 - 28] = .4;
	fjac[(j * 5 + 1) * 5 - 27] = -(UA[3 + offset] * UA[4 + offset]) / (UA[1 + offset]
                    * UA[1 + offset]);
	fjac[(j * 5 + 2) * 5 - 27] = 0.;
	fjac[(j * 5 + 3) * 5 - 27] = UA[4 + offset] / UA[1 + offset];
	fjac[(j * 5 + 4) * 5 - 27] = UA[3 + offset] / UA[1 + offset];
	fjac[(j * 5 + 5) * 5 - 27] = 0.;
	
	fjac[(j * 5 + 1) * 5 - 26] = ((UA[2 + offset] * UA[2 + offset] + UA[3 + offset] *
                    UA[3 + offset] + 
                    UA[4 + offset] * 
                    UA[4 + offset]) * .4 / (
                    UA[1 + offset] * 
                    UA[1 + offset]) - 
                    UA[5 + offset] / 
                    UA[1 + offset] * 1.4) * (
                    UA[3 + offset] / 
                    UA[1 + offset]);
	fjac[(j * 5 + 2) * 5 - 26] = UA[2 + offset] * UA[3 + offset] * -.4 / (UA[1 + offset] * UA[1 + offset]);
	fjac[(j * 5 + 3) * 5 - 26] = UA[5 + offset] / UA[1 + offset] * 1.4 - (UA[2 + offset] * UA[2 + offset] + 
                    UA[3 + offset] * 3. * 
                    UA[3 + offset] + 
                    UA[4 + offset] * 
                    UA[4 + offset]) / (
                    UA[1 + offset] * 
                    UA[1 + offset]) * 
                    .20000000000000001;
	fjac[(j * 5 + 4) * 5 - 26] = UA[3 + offset] * UA[4 + offset] * -.4 / (UA[1 + offset] * UA[1 + offset]);
	fjac[(j * 5 + 5) * 5 - 26] = UA[3 + offset] / UA[1 + offset] * 1.4;
	
	tmp1 = 1. / UA[1 + offset];
	tmp2 = tmp1 * tmp1;
	tmp3 = tmp1 * tmp2;
	
	qjac[(j * 5 + 1) * 5 - 30] = 0.;
	qjac[(j * 5 + 2) * 5 - 30] = 0.;
	qjac[(j * 5 + 3) * 5 - 30] = 0.;
	qjac[(j * 5 + 4) * 5 - 30] = 0.;
	qjac[(j * 5 + 5) * 5 - 30] = 0.;
	
	qjac[(j * 5 + 1) * 5 - 29] = -c34 * tmp2 * UA[2 + offset];
	qjac[(j * 5 + 2) * 5 - 29] = c34 * tmp1;
	qjac[(j * 5 + 3) * 5 - 29] = 0.;
	qjac[(j * 5 + 4) * 5 - 29] = 0.;
	qjac[(j * 5 + 5) * 5 - 29] = 0.;
	
	qjac[(j * 5 + 1) * 5 - 28] = -r43 * c34 * tmp2 * UA[3 + offset];
	qjac[(j * 5 + 2) * 5 - 28] = 0.;
	qjac[(j * 5 + 3) * 5 - 28] = r43 * c34 * tmp1;
	qjac[(j * 5 + 4) * 5 - 28] = 0.;
	qjac[(j * 5 + 5) * 5 - 28] = 0.;
	
	qjac[(j * 5 + 1) * 5 - 27] = -c34 * tmp2 * UA[4 + offset];
	qjac[(j * 5 + 2) * 5 - 27] = 0.;
	qjac[(j * 5 + 3) * 5 - 27] = 0.;
	qjac[(j * 5 + 4) * 5 - 27] = c34 * tmp1;
	qjac[(j * 5 + 5) * 5 - 27] = 0.;
	
	/* Computing 2nd power */
	d_1 = UA[2 + offset];
	/* Computing 2nd power */
	d_2 = UA[3 + offset];
	/* Computing 2nd power */
	d_3 = UA[4 + offset];
	qjac[(j * 5 + 1) * 5 - 26] = -(c34 - c1345) * tmp3 * (d_1 * 
                    d_1) - (r43 * c34 - c1345) * tmp3 * (d_2 * d_2) - (
                    c34 - c1345) * tmp3 * (d_3 * d_3) - c1345 * tmp2 * 
                    UA[5 + offset];
	
	qjac[(j * 5 + 2) * 5 - 26] = (c34 - c1345) * tmp2 * UA[2 + offset];
	qjac[(j * 5 + 3) * 5 - 26] = (r43 * c34 - c1345) * tmp2 * 
                    UA[3 + offset];
	qjac[(j * 5 + 4) * 5 - 26] = (c34 - c1345) * tmp2 * UA[4 + offset];
	qjac[(j * 5 + 5) * 5 - 26] = c1345 * tmp1;
      }
      
      for (l = 1; l <= 5; ++l) {
	for (m = 1; m <= 5; ++m) {
	  AA(m,l,(5*(x+Y_A[0]+z))) = 0.;
	  BA(m,l,(5*(x+Y_A[0]+z))) = 0.;
	  CA(m,l,(5*(x+Y_A[0]+z))) = 0.;
	}
	BA(l,l,(5*(x+Y_A[0]+z))) = 1.;
      }
      
      for (j = 2; j <= ny - 1; ++j) {
	YO(y, j)
        offset = 5 * (x + y + z);
	tmp1 = dt * ty1;
	tmp2 = dt * ty2;
	
	AA(1,1,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 1) * 5 - 30] - tmp1 * 
                    qjac[((j - 1) * 5 + 1) * 5 - 30] - tmp1 * dy1;
	AA(1,2,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 2) * 5 - 30] - tmp1 * 
                    qjac[((j - 1) * 5 + 2) * 5 - 30];
	AA(1,3,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 3) * 5 - 30] - tmp1 * 
                    qjac[((j - 1) * 5 + 3) * 5 - 30];
	AA(1,4,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 4) * 5 - 30] - tmp1 * 
                    qjac[((j - 1) * 5 + 4) * 5 - 30];
	AA(1,5,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 5) * 5 - 30] - tmp1 * 
                    qjac[((j - 1) * 5 + 5) * 5 - 30];
	
	AA(2,1,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 1) * 5 - 29] - tmp1 * 
                    qjac[((j - 1) * 5 + 1) * 5 - 29];
	AA(2,2,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 2) * 5 - 29] - tmp1 * 
                    qjac[((j - 1) * 5 + 2) * 5 - 29] - tmp1 * dy2;
	AA(2,3,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 3) * 5 - 29] - tmp1 * 
                    qjac[((j - 1) * 5 + 3) * 5 - 29];
	AA(2,4,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 4) * 5 - 29] - tmp1 * 
                    qjac[((j - 1) * 5 + 4) * 5 - 29];
	AA(2,5,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 5) * 5 - 29] - tmp1 * 
                    qjac[((j - 1) * 5 + 5) * 5 - 29];
	
	AA(3,1,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 1) * 5 - 28] - tmp1 * 
                    qjac[((j - 1) * 5 + 1) * 5 - 28];
	AA(3,2,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 2) * 5 - 28] - tmp1 * 
                    qjac[((j - 1) * 5 + 2) * 5 - 28];
	AA(3,3,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 3) * 5 - 28] - tmp1 * 
                    qjac[((j - 1) * 5 + 3) * 5 - 28] - tmp1 * dy3;
	AA(3,4,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 4) * 5 - 28] - tmp1 * 
                    qjac[((j - 1) * 5 + 4) * 5 - 28];
	AA(3,5,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 5) * 5 - 28] - tmp1 * 
                    qjac[((j - 1) * 5 + 5) * 5 - 28];
	
	AA(4,1,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 1) * 5 - 27] - tmp1 * 
                    qjac[((j - 1) * 5 + 1) * 5 - 27];
	AA(4,2,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 2) * 5 - 27] - tmp1 * 
                    qjac[((j - 1) * 5 + 2) * 5 - 27];
	AA(4,3,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 3) * 5 - 27] - tmp1 * 
                    qjac[((j - 1) * 5 + 3) * 5 - 27];
	AA(4,4,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 4) * 5 - 27] - tmp1 * 
                    qjac[((j - 1) * 5 + 4) * 5 - 27] - tmp1 * dy4;
	AA(4,5,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 5) * 5 - 27] - tmp1 * 
                    qjac[((j - 1) * 5 + 5) * 5 - 27];
	
	AA(5,1,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 1) * 5 - 26] - tmp1 * 
                    qjac[((j - 1) * 5 + 1) * 5 - 26];
	AA(5,2,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 2) * 5 - 26] - tmp1 * 
                    qjac[((j - 1) * 5 + 2) * 5 - 26];
	AA(5,3,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 3) * 5 - 26] - tmp1 * 
                    qjac[((j - 1) * 5 + 3) * 5 - 26];
	AA(5,4,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 4) * 5 - 26] - tmp1 * 
                    qjac[((j - 1) * 5 + 4) * 5 - 26];
	AA(5,5,offset) = 
                    -tmp2 * fjac[((j - 1) * 5 + 5) * 5 - 26] - tmp1 * 
                    qjac[((j - 1) * 5 + 5) * 5 - 26] - tmp1 * dy5;
	
	BA(1,1,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 1) * 5 - 30] + 1. + tmp1 * 2. * 
                    dy1;
	BA(1,2,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 2) * 5 - 30];
	BA(1,3,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 3) * 5 - 30];
	BA(1,4,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 4) * 5 - 30];
	BA(1,5,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 5) * 5 - 30];
	
	BA(2,1,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 1) * 5 - 29];
	BA(2,2,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 2) * 5 - 29] + 1. + tmp1 * 2. * 
                    dy2;
	BA(2,3,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 3) * 5 - 29];
	BA(2,4,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 4) * 5 - 29];
	BA(2,5,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 5) * 5 - 29];
	
	BA(3,1,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 1) * 5 - 28];
	BA(3,2,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 2) * 5 - 28];
	BA(3,3,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 3) * 5 - 28] + 1. + tmp1 * 2. * 
                    dy3;
	BA(3,4,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 4) * 5 - 28];
	BA(3,5,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 5) * 5 - 28];
	
	BA(4,1,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 1) * 5 - 27];
	BA(4,2,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 2) * 5 - 27];
	BA(4,3,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 3) * 5 - 27];
	BA(4,4,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 4) * 5 - 27] + 1. + tmp1 * 2. * 
                    dy4;
	BA(4,5,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 5) * 5 - 27];
	
	BA(5,1,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 1) * 5 - 26];
	BA(5,2,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 2) * 5 - 26];
	BA(5,3,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 3) * 5 - 26];
	BA(5,4,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 4) * 5 - 26];
	BA(5,5,offset) = tmp1 
                    * 2. * qjac[(j * 5 + 5) * 5 - 26] + 1. + tmp1 * 2. * 
                    dy5;
	
	CA(1,1,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 1) * 5 - 30] - tmp1 * qjac[((j 
                    + 1) * 5 + 1) * 5 - 30] - tmp1 * dy1;
	CA(1,2,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 2) * 5 - 30] - tmp1 * qjac[((j 
                    + 1) * 5 + 2) * 5 - 30];
	CA(1,3,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 3) * 5 - 30] - tmp1 * qjac[((j 
                    + 1) * 5 + 3) * 5 - 30];
	CA(1,4,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 4) * 5 - 30] - tmp1 * qjac[((j 
                    + 1) * 5 + 4) * 5 - 30];
	CA(1,5,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 5) * 5 - 30] - tmp1 * qjac[((j 
                    + 1) * 5 + 5) * 5 - 30];
	
	CA(2,1,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 1) * 5 - 29] - tmp1 * qjac[((j 
                    + 1) * 5 + 1) * 5 - 29];
	CA(2,2,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 2) * 5 - 29] - tmp1 * qjac[((j 
                    + 1) * 5 + 2) * 5 - 29] - tmp1 * dy2;
	CA(2,3,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 3) * 5 - 29] - tmp1 * qjac[((j 
                    + 1) * 5 + 3) * 5 - 29];
	CA(2,4,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 4) * 5 - 29] - tmp1 * qjac[((j 
                    + 1) * 5 + 4) * 5 - 29];
	CA(2,5,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 5) * 5 - 29] - tmp1 * qjac[((j 
                    + 1) * 5 + 5) * 5 - 29];
	
	CA(3,1,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 1) * 5 - 28] - tmp1 * qjac[((j 
                    + 1) * 5 + 1) * 5 - 28];
	CA(3,2,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 2) * 5 - 28] - tmp1 * qjac[((j 
                    + 1) * 5 + 2) * 5 - 28];
	CA(3,3,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 3) * 5 - 28] - tmp1 * qjac[((j 
                    + 1) * 5 + 3) * 5 - 28] - tmp1 * dy3;
	CA(3,4,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 4) * 5 - 28] - tmp1 * qjac[((j 
                    + 1) * 5 + 4) * 5 - 28];
	CA(3,5,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 5) * 5 - 28] - tmp1 * qjac[((j 
                    + 1) * 5 + 5) * 5 - 28];
	
	CA(4,1,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 1) * 5 - 27] - tmp1 * qjac[((j 
                    + 1) * 5 + 1) * 5 - 27];
	CA(4,2,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 2) * 5 - 27] - tmp1 * qjac[((j 
                    + 1) * 5 + 2) * 5 - 27];
	CA(4,3,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 3) * 5 - 27] - tmp1 * qjac[((j 
                    + 1) * 5 + 3) * 5 - 27];
	CA(4,4,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 4) * 5 - 27] - tmp1 * qjac[((j 
                    + 1) * 5 + 4) * 5 - 27] - tmp1 * dy4;
	CA(4,5,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 5) * 5 - 27] - tmp1 * qjac[((j 
                    + 1) * 5 + 5) * 5 - 27];
	
	CA(5,1,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 1) * 5 - 26] - tmp1 * qjac[((j 
                    + 1) * 5 + 1) * 5 - 26];
	CA(5,2,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 2) * 5 - 26] - tmp1 * qjac[((j 
                    + 1) * 5 + 2) * 5 - 26];
	CA(5,3,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 3) * 5 - 26] - tmp1 * qjac[((j 
                    + 1) * 5 + 3) * 5 - 26];
	CA(5,4,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 4) * 5 - 26] - tmp1 * qjac[((j 
                    + 1) * 5 + 4) * 5 - 26];
	CA(5,5,offset) = tmp2 
                    * fjac[((j + 1) * 5 + 5) * 5 - 26] - tmp1 * qjac[((j 
                    + 1) * 5 + 5) * 5 - 26] - tmp1 * dy5;
      }
      
      for (l = 1; l <= 5; ++l) {
	for (m = 1; m <= 5; ++m) {
	  AA(m,l,(5*(x+Y_A[ny-1]+z))) = 0.;
	  BA(m,l,(5*(x+Y_A[ny-1]+z))) = 0.;
	  CA(m,l,(5*(x+Y_A[ny-1]+z))) = 0.;
	}
	BA(l,l,(5*(x+Y_A[ny-1]+z))) = 1.;
      }
    }
  }
  
  return(0);
} /* jacy_ */


int jacz_(dt, tz1, tz2, dz1, dz2, dz3, dz4, dz5, nx, ny, nz)
double dt, tz1, tz2, dz1, dz2, dz3, dz4, dz5;
int nx, ny, nz;
{
/* System generated locals */
  double d_1, d_2, d_3;
  
/* Local variables */
  double fjac[5*5*S3];
  double sjac[5*5*S3];
  int i, j, k, l, m;
  double c34, r43, c1345, tmp1, tmp2, tmp3;
  int x, y, z, offset;
  
/* ***form the block tridiagonal system for zeta-direction sweep */
  
  r43 = 1.3333333333333333;
  c34 = .10000000000000001;
  c1345 = .19599999999999998;
  
  for (j = 2; j <= ny - 1; ++j) {
    YO(y, j)
    for (i = 2; i <= nx - 1; ++i) {
      XO(x, i)
/* ***compute the zeta-direction flux jacobian */
      
      for (k = 1; k <= nz; ++k) {
	ZO(z, k)
        offset = 5 * (x + y + z) - 1;
	fjac[(k * 5 + 1) * 5 - 30] = 0.;
	fjac[(k * 5 + 2) * 5 - 30] = 0.;
	fjac[(k * 5 + 3) * 5 - 30] = 0.;
	fjac[(k * 5 + 4) * 5 - 30] = 1.;
	fjac[(k * 5 + 5) * 5 - 30] = 0.;
	
	fjac[(k * 5 + 1) * 5 - 29] = -(UA[2 + offset] * UA[4 + offset]) / (UA[1 + offset]
                    * UA[1 + offset]);
	fjac[(k * 5 + 2) * 5 - 29] = UA[4 + offset] / UA[1 + offset];
	fjac[(k * 5 + 3) * 5 - 29] = 0.;
	fjac[(k * 5 + 4) * 5 - 29] = UA[2 + offset] / UA[1 + offset];
	fjac[(k * 5 + 5) * 5 - 29] = 0.;
	
	fjac[(k * 5 + 1) * 5 - 28] = -(UA[3 + offset] * UA[4 + offset]) / (UA[1 + offset]
                    * UA[1 + offset]);
	fjac[(k * 5 + 2) * 5 - 28] = 0.;
	fjac[(k * 5 + 3) * 5 - 28] = UA[4 + offset] / UA[1 + offset];
	fjac[(k * 5 + 4) * 5 - 28] = UA[3 + offset] / UA[1 + offset];
	fjac[(k * 5 + 5) * 5 - 28] = 0.;
	
	/* Computing 2nd power */
	d_1 = UA[4 + offset] / UA[1 + offset];
	fjac[(k * 5 + 1) * 5 - 27] = -(d_1 * d_1) + (UA[2 + offset] * UA[2 + offset] + UA[3 + offset] * UA[3 + offset] + UA[4 + offset] *
                    UA[4 + offset]) / (
                    UA[1 + offset] * 
                    UA[1 + offset]) * 
                    .20000000000000001;
	fjac[(k * 5 + 2) * 5 - 27] = UA[2 + offset] / UA[1 + offset] * -.4;
	fjac[(k * 5 + 3) * 5 - 27] = UA[3 + offset] / UA[1 + offset] * -.4;
	fjac[(k * 5 + 4) * 5 - 27] = UA[4 + offset] / UA[1 + offset] * 1.6000000000000001;
	fjac[(k * 5 + 5) * 5 - 27] = .4;
	
	fjac[(k * 5 + 1) * 5 - 26] = ((UA[2 + offset] * UA[2 + offset] + UA[3 + offset] *
                    UA[3 + offset] + 
                    UA[4 + offset] * 
                    UA[4 + offset]) * .4 / (
                    UA[1 + offset] * 
                    UA[1 + offset]) - 
                    UA[5 + offset] / 
                    UA[1 + offset] * 1.4) * (
                    UA[4 + offset] / 
                    UA[1 + offset]);
	fjac[(k * 5 + 2) * 5 - 26] = UA[2 + offset] * UA[4 + offset] * -.4 / (UA[1 + offset] * UA[1 + offset]);
	fjac[(k * 5 + 3) * 5 - 26] = UA[3 + offset] * UA[4 + offset] * -.4 / (UA[1 + offset] * UA[1 + offset]);
	fjac[(k * 5 + 4) * 5 - 26] = UA[5 + offset] / UA[1 + offset] * 1.4 - (UA[2 + offset] * UA[2 + offset] + 
                    UA[3 + offset] * 
                    UA[3 + offset] + 
                    UA[4 + offset] * 3. * 
                    UA[4 + offset]) / (
                    UA[1 + offset] * 
                    UA[1 + offset]) * 
                    .20000000000000001;
	fjac[(k * 5 + 5) * 5 - 26] = UA[4 + offset] / UA[1 + offset] * 1.4;
	
	tmp1 = 1. / UA[1 + offset];
	tmp2 = tmp1 * tmp1;
	tmp3 = tmp1 * tmp2;
	
	sjac[(k * 5 + 1) * 5 - 30] = 0.;
	sjac[(k * 5 + 2) * 5 - 30] = 0.;
	sjac[(k * 5 + 3) * 5 - 30] = 0.;
	sjac[(k * 5 + 4) * 5 - 30] = 0.;
	sjac[(k * 5 + 5) * 5 - 30] = 0.;
	
	sjac[(k * 5 + 1) * 5 - 29] = -c34 * tmp2 * UA[2 + offset];
	sjac[(k * 5 + 2) * 5 - 29] = c34 * tmp1;
	sjac[(k * 5 + 3) * 5 - 29] = 0.;
	sjac[(k * 5 + 4) * 5 - 29] = 0.;
	sjac[(k * 5 + 5) * 5 - 29] = 0.;
	
	sjac[(k * 5 + 1) * 5 - 28] = -c34 * tmp2 * UA[3 + offset];
	sjac[(k * 5 + 2) * 5 - 28] = 0.;
	sjac[(k * 5 + 3) * 5 - 28] = c34 * tmp1;
	sjac[(k * 5 + 4) * 5 - 28] = 0.;
	sjac[(k * 5 + 5) * 5 - 28] = 0.;
	
	sjac[(k * 5 + 1) * 5 - 27] = -r43 * c34 * tmp2 * UA[4 + offset];
	sjac[(k * 5 + 2) * 5 - 27] = 0.;
	sjac[(k * 5 + 3) * 5 - 27] = 0.;
	sjac[(k * 5 + 4) * 5 - 27] = r43 * .1 * 1. * tmp1;
	sjac[(k * 5 + 5) * 5 - 27] = 0.;
	
	/* Computing 2nd power */
	d_1 = UA[2 + offset];
	/* Computing 2nd power */
	d_2 = UA[3 + offset];
	/* Computing 2nd power */
	d_3 = UA[4 + offset];
	sjac[(k * 5 + 1) * 5 - 26] = -(c34 - c1345) * tmp3 * (d_1 * 
                    d_1) - (c34 - c1345) * tmp3 * (d_2 * d_2) - (r43 * 
                    c34 - c1345) * tmp3 * (d_3 * d_3) - c1345 * tmp2 * 
                    UA[5 + offset];
	
	sjac[(k * 5 + 2) * 5 - 26] = (c34 - c1345) * tmp2 * UA[2 + offset];
	sjac[(k * 5 + 3) * 5 - 26] = (c34 - c1345) * tmp2 * UA[3 + offset];
	sjac[(k * 5 + 4) * 5 - 26] = (r43 * c34 - c1345) * tmp2 * 
                    UA[4 + offset];
	sjac[(k * 5 + 5) * 5 - 26] = c1345 * tmp1;
      }
      
      for (l = 1; l <= 5; ++l) {
	for (m = 1; m <= 5; ++m) {
	  AA(m,l,(5*(x+y+Z_A[0]))) = 0.;
	  BA(m,l,(5*(x+y+Z_A[0]))) = 0.;
	  CA(m,l,(5*(x+y+Z_A[0]))) = 0.;
	}
	BA(l,l,(5*(x+y+Z_A[0]))) = 1.;
      }
      
      for (k = 2; k <= nz - 1; ++k) {
	ZO(z, k)
        offset = 5 * (x + y + z);
	tmp1 = dt * tz1;
	tmp2 = dt * tz2;
	
	AA(1,1,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 1) * 5 - 30] - tmp1 * 
                    sjac[((k - 1) * 5 + 1) * 5 - 30] - tmp1 * dz1;
	AA(1,2,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 2) * 5 - 30] - tmp1 * 
                    sjac[((k - 1) * 5 + 2) * 5 - 30];
	AA(1,3,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 3) * 5 - 30] - tmp1 * 
                    sjac[((k - 1) * 5 + 3) * 5 - 30];
	AA(1,4,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 4) * 5 - 30] - tmp1 * 
                    sjac[((k - 1) * 5 + 4) * 5 - 30];
	AA(1,5,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 5) * 5 - 30] - tmp1 * 
                    sjac[((k - 1) * 5 + 5) * 5 - 30];
	
	AA(2,1,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 1) * 5 - 29] - tmp1 * 
                    sjac[((k - 1) * 5 + 1) * 5 - 29];
	AA(2,2,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 2) * 5 - 29] - tmp1 * 
                    sjac[((k - 1) * 5 + 2) * 5 - 29] - tmp1 * dz2;
	AA(2,3,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 3) * 5 - 29] - tmp1 * 
                    sjac[((k - 1) * 5 + 3) * 5 - 29];
	AA(2,4,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 4) * 5 - 29] - tmp1 * 
                    sjac[((k - 1) * 5 + 4) * 5 - 29];
	AA(2,5,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 5) * 5 - 29] - tmp1 * 
                    sjac[((k - 1) * 5 + 5) * 5 - 29];
	
	AA(3,1,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 1) * 5 - 28] - tmp1 * 
                    sjac[((k - 1) * 5 + 1) * 5 - 28];
	AA(3,2,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 2) * 5 - 28] - tmp1 * 
                    sjac[((k - 1) * 5 + 2) * 5 - 28];
	AA(3,3,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 3) * 5 - 28] - tmp1 * 
                    sjac[((k - 1) * 5 + 3) * 5 - 28] - tmp1 * dz3;
	AA(3,4,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 4) * 5 - 28] - tmp1 * 
                    sjac[((k - 1) * 5 + 4) * 5 - 28];
	AA(3,5,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 5) * 5 - 28] - tmp1 * 
                    sjac[((k - 1) * 5 + 5) * 5 - 28];
	
	AA(4,1,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 1) * 5 - 27] - tmp1 * 
                    sjac[((k - 1) * 5 + 1) * 5 - 27];
	AA(4,2,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 2) * 5 - 27] - tmp1 * 
                    sjac[((k - 1) * 5 + 2) * 5 - 27];
	AA(4,3,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 3) * 5 - 27] - tmp1 * 
                    sjac[((k - 1) * 5 + 3) * 5 - 27];
	AA(4,4,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 4) * 5 - 27] - tmp1 * 
                    sjac[((k - 1) * 5 + 4) * 5 - 27] - tmp1 * dz4;
	AA(4,5,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 5) * 5 - 27] - tmp1 * 
                    sjac[((k - 1) * 5 + 5) * 5 - 27];
	
	AA(5,1,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 1) * 5 - 26] - tmp1 * 
                    sjac[((k - 1) * 5 + 1) * 5 - 26];
	AA(5,2,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 2) * 5 - 26] - tmp1 * 
                    sjac[((k - 1) * 5 + 2) * 5 - 26];
	AA(5,3,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 3) * 5 - 26] - tmp1 * 
                    sjac[((k - 1) * 5 + 3) * 5 - 26];
	AA(5,4,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 4) * 5 - 26] - tmp1 * 
                    sjac[((k - 1) * 5 + 4) * 5 - 26];
	AA(5,5,offset) = 
                    -tmp2 * fjac[((k - 1) * 5 + 5) * 5 - 26] - tmp1 * 
                    sjac[((k - 1) * 5 + 5) * 5 - 26] - tmp1 * dz5;
	
	BA(1,1,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 1) * 5 - 30] + 1. + tmp1 * 2. * 
                    dz1;
	BA(1,2,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 2) * 5 - 30];
	BA(1,3,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 3) * 5 - 30];
	BA(1,4,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 4) * 5 - 30];
	BA(1,5,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 5) * 5 - 30];
	
	BA(2,1,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 1) * 5 - 29];
	BA(2,2,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 2) * 5 - 29] + 1. + tmp1 * 2. * 
                    dz2;
	BA(2,3,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 3) * 5 - 29];
	BA(2,4,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 4) * 5 - 29];
	BA(2,5,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 5) * 5 - 29];
	
	BA(3,1,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 1) * 5 - 28];
	BA(3,2,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 2) * 5 - 28];
	BA(3,3,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 3) * 5 - 28] + 1. + tmp1 * 2. * 
                    dz3;
	BA(3,4,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 4) * 5 - 28];
	BA(3,5,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 5) * 5 - 28];
	
	BA(4,1,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 1) * 5 - 27];
	BA(4,2,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 2) * 5 - 27];
	BA(4,3,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 3) * 5 - 27];
	BA(4,4,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 4) * 5 - 27] + 1. + tmp1 * 2. * 
                    dz4;
	BA(4,5,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 5) * 5 - 27];
	
	BA(5,1,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 1) * 5 - 26];
	BA(5,2,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 2) * 5 - 26];
	BA(5,3,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 3) * 5 - 26];
	BA(5,4,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 4) * 5 - 26];
	BA(5,5,offset) = tmp1 
                    * 2. * sjac[(k * 5 + 5) * 5 - 26] + 1. + tmp1 * 2. * 
                    dz5;
	
	CA(1,1,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 1) * 5 - 30] - tmp1 * sjac[((k 
                    + 1) * 5 + 1) * 5 - 30] - tmp1 * dz1;
	CA(1,2,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 2) * 5 - 30] - tmp1 * sjac[((k 
                    + 1) * 5 + 2) * 5 - 30];
	CA(1,3,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 3) * 5 - 30] - tmp1 * sjac[((k 
                    + 1) * 5 + 3) * 5 - 30];
	CA(1,4,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 4) * 5 - 30] - tmp1 * sjac[((k 
                    + 1) * 5 + 4) * 5 - 30];
	CA(1,5,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 5) * 5 - 30] - tmp1 * sjac[((k 
                    + 1) * 5 + 5) * 5 - 30];
	
	CA(2,1,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 1) * 5 - 29] - tmp1 * sjac[((k 
                    + 1) * 5 + 1) * 5 - 29];
	CA(2,2,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 2) * 5 - 29] - tmp1 * sjac[((k 
                    + 1) * 5 + 2) * 5 - 29] - tmp1 * dz2;
	CA(2,3,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 3) * 5 - 29] - tmp1 * sjac[((k 
                    + 1) * 5 + 3) * 5 - 29];
	CA(2,4,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 4) * 5 - 29] - tmp1 * sjac[((k 
                    + 1) * 5 + 4) * 5 - 29];
	CA(2,5,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 5) * 5 - 29] - tmp1 * sjac[((k 
                    + 1) * 5 + 5) * 5 - 29];
	
	CA(3,1,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 1) * 5 - 28] - tmp1 * sjac[((k 
                    + 1) * 5 + 1) * 5 - 28];
	CA(3,2,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 2) * 5 - 28] - tmp1 * sjac[((k 
                    + 1) * 5 + 2) * 5 - 28];
	CA(3,3,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 3) * 5 - 28] - tmp1 * sjac[((k 
                    + 1) * 5 + 3) * 5 - 28] - tmp1 * dz3;
	CA(3,4,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 4) * 5 - 28] - tmp1 * sjac[((k 
                    + 1) * 5 + 4) * 5 - 28];
	CA(3,5,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 5) * 5 - 28] - tmp1 * sjac[((k 
                    + 1) * 5 + 5) * 5 - 28];
	
	CA(4,1,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 1) * 5 - 27] - tmp1 * sjac[((k 
                    + 1) * 5 + 1) * 5 - 27];
	CA(4,2,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 2) * 5 - 27] - tmp1 * sjac[((k 
                    + 1) * 5 + 2) * 5 - 27];
	CA(4,3,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 3) * 5 - 27] - tmp1 * sjac[((k 
                    + 1) * 5 + 3) * 5 - 27];
	CA(4,4,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 4) * 5 - 27] - tmp1 * sjac[((k 
                    + 1) * 5 + 4) * 5 - 27] - tmp1 * dz4;
	CA(4,5,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 5) * 5 - 27] - tmp1 * sjac[((k 
                    + 1) * 5 + 5) * 5 - 27];
	
	CA(5,1,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 1) * 5 - 26] - tmp1 * sjac[((k 
                    + 1) * 5 + 1) * 5 - 26];
	CA(5,2,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 2) * 5 - 26] - tmp1 * sjac[((k 
                    + 1) * 5 + 2) * 5 - 26];
	CA(5,3,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 3) * 5 - 26] - tmp1 * sjac[((k 
                    + 1) * 5 + 3) * 5 - 26];
	CA(5,4,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 4) * 5 - 26] - tmp1 * sjac[((k 
                    + 1) * 5 + 4) * 5 - 26];
	CA(5,5,offset) = tmp2 
                    * fjac[((k + 1) * 5 + 5) * 5 - 26] - tmp1 * sjac[((k 
                    + 1) * 5 + 5) * 5 - 26] - tmp1 * dz5;
      }
      
      for (l = 1; l <= 5; ++l) {
	for (m = 1; m <= 5; ++m) {
	  AA(m,l,(5*(x+y+Z_A[nz-1]))) = 0.;
	  BA(m,l,(5*(x+y+Z_A[nz-1]))) = 0.;
	  CA(m,l,(5*(x+y+Z_A[nz-1]))) = 0.;
	}
	BA(l,l,(5*(x+y+Z_A[nz-1]))) = 1.;
      }
    }
  }
  
  return(0);
} /* jacz_ */
