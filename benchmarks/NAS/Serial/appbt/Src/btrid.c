/* --- btrid.c ---
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

int btridx_(nx, ny, nz)
int nx, ny, nz;
{
  int i, j, k, m, x, y, z, offset;
  double tm[25], tv[5], tmp;
  
/****solution of multiple, independent systems of block tridiagonal system
    s*/
/*   using Gaussian elimination (without pivoting) algorithm */
/*   (block size = 5 X 5) */
  
/* ***forward elimination phase */
  
  for (k = 2; k <= nz - 1; ++k) {
    ZO(z, k)
    for (j = 2; j <= ny - 1; ++j) {
      YO(y, j)
      offset = 5 * (X_A[0] + y + z);
      tmp = 1. / BA(1,1,offset);
      
      BA(2,1,offset) *= tmp;
      BA(3,1,offset) *= tmp;
      BA(4,1,offset) *= tmp;
      BA(5,1,offset) *= tmp;
      
      BA(2,2,offset) -= BA(2,1,offset) * BA(1,2,offset);
      
      tmp = 1. / BA(2,2,offset);
      
      BA(3,2,offset) = (BA(3,2,offset) - BA(3,1,offset) * BA(1,2,offset)) * tmp;
      BA(4,2,offset) = (BA(4,2,offset) - BA(4,1,offset) * BA(1,2,offset)) * tmp;
      BA(5,2,offset) = (BA(5,2,offset) - BA(5,1,offset) * BA(1,2,offset)) * tmp;
      
      BA(2,3,offset) -= BA(2,1,offset) * BA(1,3,offset);
      BA(3,3,offset) = BA(3,3,offset) - BA(3,1,offset) * BA(1,3,offset) 
	                              - BA(3,2,offset) * BA(2,3,offset);
      
      tmp = 1. / BA(3,3,offset);
      
      BA(4,3,offset) = (BA(4,3,offset) - BA(4,1,offset) * BA(1,3,offset) - BA(4,2,offset) * BA(2,3,offset)) * tmp;
      BA(5,3,offset) = (BA(5,3,offset) - BA(5,1,offset) * BA(1,3,offset) - BA(5,2,offset) * BA(2,3,offset)) * tmp;
      
      BA(2,4,offset) -= BA(2,1,offset) * BA(1,4,offset);
      BA(3,4,offset) = BA(3,4,offset) - BA(3,1,offset) * BA(1,4,offset) - BA(3,2,offset) * BA(2,4,offset);
      BA(4,4,offset) = BA(4,4,offset) - BA(4,1,offset) * BA(1,4,offset) - BA(4,2,offset) * BA(2,4,offset) - BA(4,3,offset) * BA(3,4,offset);
      
      BA(5,4,offset) = (BA(5,4,offset) - BA(5,1,offset) * BA(1,4,offset) - BA(5,2,offset) * BA(2,4,offset) - BA(5,3,offset) * BA(3,4,offset)) / BA(4,4,offset);
      
      BA(2,5,offset) -= BA(2,1,offset) * BA(1,5,offset);
      BA(3,5,offset) = BA(3,5,offset) - BA(3,1,offset) * BA(1,5,offset) - BA(3,2,offset) * BA(2,5,offset);
      BA(4,5,offset) = BA(4,5,offset) - BA(4,1,offset) * BA(1,5,offset) - BA(4,2,offset) * BA(2,5,offset) - BA(4,3,offset) * BA(3,5,offset);
      
      BA(5,5,offset) = BA(5,5,offset) - BA(5,1,offset) * BA(1,5,offset) - BA(5,2,offset) * BA(2,5,offset) - BA(5,3,offset) * BA(3,5,offset) - BA(5,4,offset) * BA(4,5,offset);
      
      for (i = 2; i <= nx; ++i) {
	offset = 5 * (X_A[i-2] + y + z);
	for (m = 1; m <= 5; ++m) {
	  tm[m * 5 - 5] = CA(1,m,offset);
	  tm[m * 5 - 4] = CA(2,m,offset) - BA(2,1,offset) * tm[m * 5 - 5];
	  tm[m * 5 - 3] = CA(3,m,offset) - BA(3,1,offset) * tm[m * 5 - 5] - BA(3,2,offset)
                    * tm[m * 5 - 4];
	  tm[m * 5 - 2] = CA(4,m,offset) - BA(4,1,offset) * tm[m * 5 - 5] - BA(4,2,offset)
                    * tm[m * 5 - 4] - BA(4,3,offset) * tm[m * 5 - 3];
	  tm[m * 5 - 1] = CA(5,m,offset) - BA(5,1,offset) * tm[m * 5 - 5] - BA(5,2,offset)
                    * tm[m * 5 - 4] - BA(5,3,offset) * tm[m * 5 - 3] - BA(5,4,offset)
                    * tm[m * 5 - 2];
	}
	
	for (m = 1; m <= 5; ++m) {
	  
	  tm[m * 5 - 1] /= BA(5,5,offset);
	  tm[m * 5 - 2] = (tm[m * 5 - 2] - BA(4,5,offset) * tm[m * 5 - 
                    1]) / BA(4,4,offset);
	  tm[m * 5 - 3] = (tm[m * 5 - 3] - BA(3,4,offset) * tm[m * 5 - 
                    2] - BA(3,5,offset) * tm[m * 5 - 1]) / BA(3,3,offset);
	  tm[m * 5 - 4] = (tm[m * 5 - 4] - BA(2,3,offset) * tm[m * 5 - 
                    3] - BA(2,4,offset) * tm[m * 5 - 2] - BA(2,5,offset) * tm[m * 5 - 
                    1]) / BA(2,2,offset);
	  tm[m * 5 - 5] = (tm[m * 5 - 5] - BA(1,2,offset) * tm[m * 5 - 
                    4] - BA(1,3,offset) * tm[m * 5 - 3] - BA(1,4,offset) * tm[m * 5 - 
                    2] - BA(1,5,offset) * tm[m * 5 - 1]) / BA(1,1,offset);
	}
	
	XO(x, i)
	offset = 5 * (x + y + z);
	for (m = 1; m <= 5; ++m) {
	  
	  BA(1,m,offset) = BA(1,m,offset) - AA(1,1,offset) * tm[m * 5 - 5] 
	            - AA(1,2,offset) * tm[m * 5 - 4] - AA(1,3,offset) *  
                    tm[m * 5 - 3] - AA(1,4,offset) * tm[m * 5 - 2] - 
		    AA(1,5,offset) * tm[m * 5 - 1];
	  
	  BA(2,m,offset) = BA(2,m,offset) - AA(2,1,offset) * tm[m * 5 - 5] - 
	            AA(2,2,offset) * tm[m * 5 - 4] - AA(2,3,offset) * tm[m * 5 - 3] 
                   - AA(2,4,offset) * tm[m * 5 - 2] - AA(2,5,offset) * tm[m * 5 - 1];
	  
	  BA(3,m,offset) = BA(3,m,offset) - AA(3,1,offset) * tm[m * 5 - 5] - 
                    AA(3,2,offset) * tm[m * 5 - 4] - AA(3,3,offset) * tm[m * 5 - 3] 
                   - AA(3,4,offset) * tm[m * 5 - 2] - AA(3,5,offset) * tm[m * 5 - 1];
	  
	  BA(4,m,offset) = BA(4,m,offset) - AA(4,1,offset) * tm[m * 5 - 5] - 
                    AA(4,2,offset) * tm[m * 5 - 4] - AA(4,3,offset) * tm[m * 5 - 3]
                   - AA(4,4,offset) * tm[m * 5 - 2] - AA(4,5,offset) * tm[m * 5 - 1];
	  
	  BA(5,m,offset) = BA(5,m,offset) - AA(5,1,offset) * tm[m * 5 - 5] - 
                    AA(5,2,offset) * tm[m * 5 - 4] - AA(5,3,offset) * tm[m * 5 - 3] -
                    AA(5,4,offset) * tm[m * 5 - 2] - AA(5,5,offset) * tm[m * 5 - 1];
	}

	offset = 5 * (X_A[i-2] + y + z);
	tv[0] = RSD[1 + offset - 1];
	tv[1] = RSD[2 + offset - 1] - BA(2,1,offset) *
                    tv[0];
	tv[2] = RSD[3 + offset - 1] - BA(3,1,offset) *
                    tv[0] - BA(3,2,offset) * tv[1];
	tv[3] = RSD[4 + offset - 1] - BA(4,1,offset) *
                    tv[0] - BA(4,2,offset) * tv[1] - BA(4,3,offset) * tv[2];
	tv[4] = RSD[5 + offset - 1] - BA(5,1,offset) *
                    tv[0] - BA(5,2,offset) * tv[1] - BA(5,3,offset) * tv[2] - BA(5,4,offset) * tv[3];
	
	tv[4] /= BA(5,5,offset);
	tv[3] = (tv[3] - BA(4,5,offset) * tv[4]) / BA(4,4,offset);
	tv[2] = (tv[2] - BA(3,4,offset) * tv[3] - BA(3,5,offset) * tv[4]) / BA(3,3,offset);
	tv[1] = (tv[1] - BA(2,3,offset) * tv[2] - BA(2,4,offset) * tv[3] - BA(2,5,offset) * tv[4]) / BA(2,2,offset);
	
	tv[0] = (tv[0] - BA(1,2,offset) * tv[1] - BA(1,3,offset) * tv[2] - BA(1,4,offset) * tv[3] - BA(1,5,offset) *
                    tv[4]) / BA(1,1,offset);
	
        offset = 5 * (x + y + z);
	for (m = 1; m <= 5; ++m)
	  RSD[m + offset - 1] = RSD[m + offset - 1] - AA(m,1,offset) * tv[0] - AA(m,2,offset) *
                    tv[1] - AA(m,3,offset) * tv[2] - AA(m,4,offset) * tv[3] - AA(m,5,offset) * tv[4];
	
	tmp = 1. / BA(1,1,offset);
	
	BA(2,1,offset) *= tmp;
	BA(3,1,offset) *= tmp;
	BA(4,1,offset) *= tmp;
	BA(5,1,offset) *= tmp;
	BA(2,2,offset) -= BA(2,1,offset) * BA(1,2,offset);
                    
	tmp = 1. / BA(2,2,offset);
	
	BA(3,2,offset) = (BA(3,2,offset) - BA(3,1,offset) * BA(1,2,offset)) * tmp;
	BA(4,2,offset) = (BA(4,2,offset) - BA(4,1,offset) * BA(1,2,offset)) * tmp;
	BA(5,2,offset) = (BA(5,2,offset) - BA(5,1,offset) * BA(1,2,offset)) * tmp;
	
	BA(2,3,offset) -= BA(2,1,offset) * BA(1,3,offset);
	BA(3,3,offset) = BA(3,3,offset) - BA(3,1,offset) * BA(1,3,offset) - 
                         BA(3,2,offset) * BA(2,3,offset);
                    
	tmp = 1. / BA(3,3,offset);
	
	BA(4,3,offset) = (BA(4,3,offset) - BA(4,1,offset) * BA(1,3,offset) - 
                         BA(4,2,offset) * BA(2,3,offset)) * tmp;
	BA(5,3,offset) = (BA(5,3,offset) - BA(5,1,offset) * BA(1,3,offset) - 
                         BA(5,2,offset) * BA(2,3,offset)) * tmp;
	
	BA(2,4,offset) -= BA(2,1,offset) * BA(1,4,offset);
	BA(3,4,offset) = BA(3,4,offset) - BA(3,1,offset) * BA(1,4,offset) - 
                         BA(3,2,offset) * BA(2,4,offset);
	BA(4,4,offset) = BA(4,4,offset) - BA(4,1,offset) * BA(1,4,offset) - 
                   BA(4,2,offset) * BA(2,4,offset) - BA(4,3,offset) * BA(3,4,offset);
	
	BA(5,4,offset) = (BA(5,4,offset) - BA(5,1,offset) * BA(1,4,offset) - 
                         BA(5,2,offset) * BA(2,4,offset) - BA(5,3,offset) 
                         * BA(3,4,offset)) / BA(4,4,offset);
	
	BA(2,5,offset) -= BA(2,1,offset) * BA(1,5,offset);
	BA(3,5,offset) = BA(3,5,offset) - BA(3,1,offset) * BA(1,5,offset) - 
                         BA(3,2,offset) * BA(2,5,offset);
	BA(4,5,offset) = BA(4,5,offset) - BA(4,1,offset) * BA(1,5,offset) - 
                  BA(4,2,offset) * BA(2,5,offset) - BA(4,3,offset) * BA(3,5,offset);
	
	BA(5,5,offset) = BA(5,5,offset) - BA(5,1,offset) * BA(1,5,offset) - 
                         BA(5,2,offset) * BA(2,5,offset) - BA(5,3,offset)
                         * BA(3,5,offset) - BA(5,4,offset) * BA(4,5,offset);
      }
    }
  }

  /* ***back-substitution phase */
  
  for (k = 2; k <= nz - 1; ++k) {
    ZO(z, k)
    for (j = 2; j <= ny - 1; ++j) {
      YO(y, j)
      offset = 5 * (X_A[nx-1] + y + z);
      RSD[2 + offset - 1] -= BA(2,1,offset) * RSD[1 + offset - 1];
      RSD[3 + offset - 1] = RSD[3 + offset - 1] - BA(3,1,offset) * RSD[1 + offset - 1] - BA(3,2,offset) * RSD[2 + offset - 1];
      RSD[4 + offset - 1] = RSD[4 + offset - 1] - BA(4,1,offset) * RSD[1 + offset - 1] - BA(4,2,offset) * RSD[2 + offset - 1] - BA(4,3,offset) * RSD[3 + offset - 1];
      RSD[5 + offset - 1] = RSD[5 + offset - 1] - BA(5,1,offset) * RSD[1 + offset - 1] - BA(5,2,offset) * RSD[2 + offset - 1] - BA(5,3,offset) * RSD[3 + offset - 1] - BA(5,4,offset) * RSD[4 + offset - 1];
      
      RSD[5 + offset - 1] /= BA(5,5,offset);
      RSD[4 + offset - 1] = (RSD[4 + offset - 1] - BA(4,5,offset) * RSD[5 + offset - 1]) / BA(4,4,offset);
      RSD[3 + offset - 1] = (RSD[3 + offset - 1] - BA(3,4,offset) * RSD[4 + offset - 1] - BA(3,5,offset) * RSD[5 + offset - 1]) / BA(3,3,offset);
      RSD[2 + offset - 1] = (RSD[2 + offset - 1] - BA(2,3,offset) * RSD[3 + offset - 1] - BA(2,4,offset) * RSD[4 + offset - 1] - BA(2,5,offset) * RSD[5 + offset - 1]) / BA(2,2,offset);
      RSD[1 + offset - 1] = (RSD[1 + offset - 1] - BA(1,2,offset) * RSD[2 + offset - 1] - BA(1,3,offset) * RSD[3 + offset - 1] - BA(1,4,offset) * RSD[4 + offset - 1] - BA(1,5,offset) * RSD[5 + offset - 1]) / BA(1,1,offset);
      
      for (i = nx - 1; i >= 1; --i) {
	XO(x, i)
        offset = 5 * (x + y + z);
	for (m = 1; m <= 5; ++m) {
	  RSD[m + offset - 1] = RSD[m + offset - 1] - CA(m,1,offset) * RSD[1 + 5 * (X_A[i] + y + z) - 1] - CA(m,2,offset) * RSD[2 + 5 * (X_A[i] + y + z) - 1] - CA(m,3,offset) * 
                    RSD[3 + 5 * (X_A[i] + y + z) - 1] - CA(m,4,offset) *
                    RSD[4 + 5 * (X_A[i] + y + z) - 1] - 
                    CA(m,5,offset) * RSD[5 + 5 * (X_A[i] + y + z) - 1];
	}
	
	RSD[2 + offset - 1] -= BA(2,1,offset) * RSD[1 + offset - 1];
	RSD[3 + offset - 1] = RSD[3 + offset - 1] - BA(3,1,offset) * RSD[1 + offset - 1] - BA(3,2,offset) * RSD[2 + offset - 1];
	RSD[4 + offset - 1] = RSD[4 + offset - 1] - BA(4,1,offset) * RSD[1 + offset - 1] - BA(4,2,offset) * RSD[2 + offset - 1] - BA(4,3,offset) * RSD[3 + offset - 1];
	RSD[5 + offset - 1] = RSD[5 + offset - 1] - BA(5,1,offset) * RSD[1 + offset - 1] - BA(5,2,offset) * RSD[2 + offset - 1] - BA(5,3,offset) * RSD[3 + offset - 1] - BA(5,4,offset) * RSD[4 + offset - 1] ;
	
	RSD[5 + offset - 1] /= BA(5,5,offset);
	RSD[4 + offset - 1] = (RSD[4 + offset - 1] - BA(4,5,offset) * RSD[5 + offset - 1]) / BA(4,4,offset);
	RSD[3 + offset - 1] = (RSD[3 + offset - 1] - BA(3,4,offset) * RSD[4 + offset - 1] - BA(3,5,offset) * RSD[5 + offset - 1]) / BA(3,3,offset);
	RSD[2 + offset - 1] = (RSD[2 + offset - 1] - BA(2,3,offset) * RSD[3 + offset - 1] - BA(2,4,offset) * RSD[4 + offset - 1] - BA(2,5,offset) * RSD[5 + offset - 1]) / BA(2,2,offset);
	RSD[1 + offset - 1] = (RSD[1 + offset - 1] - BA(1,2,offset) * RSD[2 + offset - 1] - BA(1,3,offset) * RSD[3 + offset - 1] - BA(1,4,offset) * RSD[4 + offset - 1] - BA(1,5,offset) * RSD[5 + offset - 1]) / BA(1,1,offset);
      }
    }
  }
  
  return(0);
} /* btridx_ */


int btridy_(nx, ny, nz)
int nx, ny, nz;
{
  int i, j, k, m, x, y, z, offset;
  double tm[25], tv[5], tmp;
  
/****solution of multiple, independent systems of block tridiagonal systems*/
/*   using Gaussian elimination (without pivoting) algorithm */
/*   (block size = 5 X 5) */
  
/* ***forward elimination phase */
  
  for (k = 2; k <= nz - 1; ++k) {
    ZO(z, k)
    for (i = 2; i <= nx - 1; ++i) {
      XO(x, i)
      offset = 5 * (x + Y_A[0] + z);
      tmp = 1. / BA(1,1,offset);
      
      BA(2,1,offset) *= tmp;
      BA(3,1,offset) *= tmp;
      BA(4,1,offset) *= tmp;
      BA(5,1,offset) *= tmp;
      
      BA(2,2,offset) -= BA(2,1,offset) * BA(1,2,offset);
      
      tmp = 1. / BA(2,2,offset);
      
      BA(3,2,offset) = (BA(3,2,offset) - BA(3,1,offset) * BA(1,2,offset)) * tmp;
      BA(4,2,offset) = (BA(4,2,offset) - BA(4,1,offset) * BA(1,2,offset)) * tmp;
      BA(5,2,offset) = (BA(5,2,offset) - BA(5,1,offset) * BA(1,2,offset)) * tmp;
      
      BA(2,3,offset) -= BA(2,1,offset) * BA(1,3,offset);
      BA(3,3,offset) = BA(3,3,offset) - BA(3,1,offset) * BA(1,3,offset) - BA(3,2,offset) * BA(2,3,offset);
      
      tmp = 1. / BA(3,3,offset);
      
      BA(4,3,offset) = (BA(4,3,offset) - BA(4,1,offset) * BA(1,3,offset) - BA(4,2,offset) * BA(2,3,offset)) * tmp;
      BA(5,3,offset) = (BA(5,3,offset) - BA(5,1,offset) * BA(1,3,offset) - BA(5,2,offset) * BA(2,3,offset)) * tmp;
      
      BA(2,4,offset) -= BA(2,1,offset) * BA(1,4,offset);
      BA(3,4,offset) = BA(3,4,offset) - BA(3,1,offset) * BA(1,4,offset) - BA(3,2,offset) * BA(2,4,offset);
      BA(4,4,offset) = BA(4,4,offset) - BA(4,1,offset) * BA(1,4,offset) - BA(4,2,offset) * BA(2,4,offset) - BA(4,3,offset) * BA(3,4,offset);
      
      BA(5,4,offset) = (BA(5,4,offset) - BA(5,1,offset) * BA(1,4,offset) - BA(5,2,offset) * BA(2,4,offset) - BA(5,3,offset) * BA(3,4,offset)) / BA(4,4,offset);
      
      BA(2,5,offset) -= BA(2,1,offset) * BA(1,5,offset);
      BA(3,5,offset) = BA(3,5,offset) - BA(3,1,offset) * BA(1,5,offset) - BA(3,2,offset) * BA(2,5,offset);
      BA(4,5,offset) = BA(4,5,offset) - BA(4,1,offset) * BA(1,5,offset) - BA(4,2,offset) * BA(2,5,offset) - BA(4,3,offset) * BA(3,5,offset);
      
      BA(5,5,offset) = BA(5,5,offset) - BA(5,1,offset) * BA(1,5,offset) - BA(5,2,offset) * BA(2,5,offset) - BA(5,3,offset) * BA(3,5,offset) - BA(5,4,offset) * BA(4,5,offset);
      
      for (j = 2; j <= ny; ++j) {
        YO(y, j-1)
        offset = 5 * (x + y + z);
	for (m = 1; m <= 5; ++m) {
	  
	  tm[m * 5 - 5] = CA(1,m,offset);
	  tm[m * 5 - 4] = CA(2,m,offset) - BA(2,1,offset) * tm[m * 5 - 5];
	  tm[m * 5 - 3] = CA(3,m,offset) - BA(3,1,offset) * tm[m * 5 - 5] - BA(3,2,offset)
                    * tm[m * 5 - 4];
	  tm[m * 5 - 2] = CA(4,m,offset) - BA(4,1,offset) * tm[m * 5 - 5] - BA(4,2,offset)
                    * tm[m * 5 - 4] - BA(4,3,offset) * tm[m * 5 - 3];
	  tm[m * 5 - 1] = CA(5,m,offset) - BA(5,1,offset) * tm[m * 5 - 5] - BA(5,2,offset)
                    * tm[m * 5 - 4] - BA(5,3,offset) * tm[m * 5 - 3] - BA(5,4,offset)
                    * tm[m * 5 - 2];
	}
	
	for (m = 1; m <= 5; ++m) {
	  
	  tm[m * 5 - 1] /= BA(5,5,offset);
	  tm[m * 5 - 2] = (tm[m * 5 - 2] - BA(4,5,offset) * tm[m * 5 - 
                    1]) / BA(4,4,offset);
	  tm[m * 5 - 3] = (tm[m * 5 - 3] - BA(3,4,offset) * tm[m * 5 - 
                    2] - BA(3,5,offset) * tm[m * 5 - 1]) / BA(3,3,offset);
	  tm[m * 5 - 4] = (tm[m * 5 - 4] - BA(2,3,offset) * tm[m * 5 - 
                    3] - BA(2,4,offset) * tm[m * 5 - 2] - BA(2,5,offset) * tm[m * 5 - 
                    1]) / BA(2,2,offset);
	  tm[m * 5 - 5] = (tm[m * 5 - 5] - BA(1,2,offset) * tm[m * 5 - 
                    4] - BA(1,3,offset) * tm[m * 5 - 3] - BA(1,4,offset) * tm[m * 5 - 
                    2] - BA(1,5,offset) * tm[m * 5 - 1]) / BA(1,1,offset);
	}
	
        YO(y, j);
        offset = 5 * (x + y + z);
	for (m = 1; m <= 5; ++m) {
	  
	  BA(1,m,offset) = BA(1,m,offset) 
                    - AA(1,1,offset) * tm[m * 5 - 5] - AA(1,2,offset) * tm[m * 5 - 4] - AA(1,3,offset) * 
                    tm[m * 5 - 3] - AA(1,4,offset) * tm[m * 5 - 2] - AA(1,5,offset) * tm[m * 5 - 1];
	  
	  BA(2,m,offset) = BA(2,m,offset) 
                    - AA(2,1,offset) * tm[m * 5 - 5] - AA(2,2,offset) * tm[m * 5 - 4] - AA(2,3,offset) * 
                    tm[m * 5 - 3] - AA(2,4,offset) * tm[m * 5 - 2] - AA(2,5,offset) * tm[m * 
                    5 - 1];
	  
	  BA(3,m,offset) = BA(3,m,offset) 
                    - AA(3,1,offset) * tm[m * 5 - 5] - AA(3,2,offset) * tm[m * 5 - 4] - AA(3,3,offset) * 
                    tm[m * 5 - 3] - AA(3,4,offset) * tm[m * 5 - 2] - AA(3,5,offset) * tm[m * 
                    5 - 1];
	  
	  BA(4,m,offset) = BA(4,m,offset) 
                    - AA(4,1,offset) * tm[m * 5 - 5] - AA(4,2,offset) * tm[m * 5 - 4] - AA(4,3,offset) * 
                    tm[m * 5 - 3] - AA(4,4,offset) * tm[m * 5 - 2] - AA(4,5,offset) * tm[m * 
                    5 - 1];
	  
	  BA(5,m,offset) = BA(5,m,offset) 
                    - AA(5,1,offset) * tm[m * 5 - 5] - AA(5,2,offset) * tm[m * 5 - 4] - AA(5,3,offset) * 
                    tm[m * 5 - 3] - AA(5,4,offset) * tm[m * 5 - 2] - AA(5,5,offset) * tm[m * 
                    5 - 1];
	}
	
        offset = 5 * (x + Y_A[j-2] + z);
	tv[0] = RSD[1 + offset - 1];
	tv[1] = RSD[2 + offset - 1] - BA(2,1,offset) *
                    tv[0];
	tv[2] = RSD[3 + offset - 1] - BA(3,1,offset) *
                    tv[0] - BA(3,2,offset) * tv[1];
	tv[3] = RSD[4 + offset - 1] - BA(4,1,offset) *
                    tv[0] - BA(4,2,offset) * tv[1] - BA(4,3,offset) * tv[2];
	tv[4] = RSD[5 + offset - 1] - BA(5,1,offset) *
                    tv[0] - BA(5,2,offset) * tv[1] - BA(5,3,offset) * tv[2] - BA(5,4,offset) * tv[3];
	
	tv[4] /= BA(5,5,offset);
	tv[3] = (tv[3] - BA(4,5,offset) * tv[4]) / BA(4,4,offset);
	tv[2] = (tv[2] - BA(3,4,offset) * tv[3] - BA(3,5,offset) * tv[4]) / BA(3,3,offset);
	tv[1] = (tv[1] - BA(2,3,offset) * tv[2] - BA(2,4,offset) * tv[3] - BA(2,5,offset) * tv[4]) / BA(2,2,offset);
	
	tv[0] = (tv[0] - BA(1,2,offset) * tv[1] - BA(1,3,offset) * tv[2] - BA(1,4,offset) * tv[3] - BA(1,5,offset) *
                    tv[4]) / BA(1,1,offset);
	
        offset = 5 * (x + y + z);
	for (m = 1; m <= 5; ++m)
	  RSD[m + offset - 1] = RSD[m + offset - 1] - AA(m,1,offset) * tv[0] - AA(m,2,offset) *
                    tv[1] - AA(m,3,offset) * tv[2] - AA(m,4,offset) * tv[3] - AA(m,5,offset) * tv[4];
	
	tmp = 1. / BA(1,1,offset);
	
	BA(2,1,offset) *= tmp;
	BA(3,1,offset) *= tmp;
	BA(4,1,offset) *= tmp;
	BA(5,1,offset) *= tmp;
	
	BA(2,2,offset) -= BA(2,1,offset) * BA(1,2,offset);
	
	tmp = 1. / BA(2,2,offset);
	
	BA(3,2,offset) = (BA(3,2,offset) - BA(3,1,offset) * BA(1,2,offset)) * tmp;
	BA(4,2,offset) = (BA(4,2,offset) - BA(4,1,offset) * BA(1,2,offset)) * tmp;
	BA(5,2,offset) = (BA(5,2,offset) - BA(5,1,offset) * BA(1,2,offset)) * tmp;
	
	BA(2,3,offset) -= BA(2,1,offset) * BA(1,3,offset);
	BA(3,3,offset) = BA(3,3,offset) - BA(3,1,offset) * BA(1,3,offset) - BA(3,2,offset) * 
                    BA(2,3,offset);
	
	tmp = 1. / BA(3,3,offset);
	
	BA(4,3,offset) = (BA(4,3,offset) - BA(4,1,offset) * BA(1,3,offset) - BA(4,2,offset) * 
                    BA(2,3,offset)) 
                    * tmp;
	BA(5,3,offset) = (BA(5,3,offset) - BA(5,1,offset) * BA(1,3,offset) - BA(5,2,offset) * 
                    BA(2,3,offset)) 
                    * tmp;
	
	BA(2,4,offset) -= BA(2,1,offset) * BA(1,4,offset);
	BA(3,4,offset) = BA(3,4,offset) - BA(3,1,offset) * BA(1,4,offset) - BA(3,2,offset) * 
                    BA(2,4,offset);
	BA(4,4,offset) = BA(4,4,offset) - BA(4,1,offset) * BA(1,4,offset) - BA(4,2,offset) * 
                    BA(2,4,offset) 
                    - BA(4,3,offset)
                    * BA(3,4,offset);
	
	BA(5,4,offset) = (BA(5,4,offset) - BA(5,1,offset) * BA(1,4,offset) - BA(5,2,offset) * 
                    BA(2,4,offset) 
                    - BA(5,3,offset)
                    * BA(3,4,offset)) / BA(4,4,offset);
	
	BA(2,5,offset) -= BA(2,1,offset) * BA(1,5,offset);
        BA(3,5,offset) = BA(3,5,offset) - BA(3,1,offset) * BA(1,5,offset) - BA(3,2,offset) * 
                    BA(2,5,offset);
	BA(4,5,offset) = BA(4,5,offset) - BA(4,1,offset) * BA(1,5,offset) - BA(4,2,offset) * 
                    BA(2,5,offset) 
                    - BA(4,3,offset)
                    * BA(3,5,offset);
	
	BA(5,5,offset) = BA(5,5,offset) - BA(5,1,offset) * BA(1,5,offset) - BA(5,2,offset) * 
                    BA(2,5,offset) 
                    - BA(5,3,offset)
                    * BA(3,5,offset) - BA(5,4,offset) * BA(4,5,offset);
      }
    }
  }
  
/* ***back-substitution phase */
  
  for (k = 2; k <= nz - 1; ++k) {
    ZO(z, k)
    for (i = 2; i <= nx - 1; ++i) {
      XO(x, i)
      offset = 5 * (x + Y_A[ny-1] + z);
      RSD[2 + offset - 1] -= BA(2,1,offset) * RSD[1 + offset - 1];
      RSD[3 + offset - 1] = RSD[3 + offset - 1] - BA(3,1,offset) * RSD[1 + offset - 1] - BA(3,2,offset) * RSD[2 + offset - 1];
      RSD[4 + offset - 1] = RSD[4 + offset - 1] - BA(4,1,offset) * RSD[1 + offset - 1] - BA(4,2,offset) * RSD[2 + offset - 1] - BA(4,3,offset) * RSD[3 + offset - 1];
      RSD[5 + offset - 1] = RSD[5 + offset - 1] - BA(5,1,offset) * RSD[1 + offset - 1] - BA(5,2,offset) * RSD[2 + offset - 1] - BA(5,3,offset) * RSD[3 + offset - 1] - BA(5,4,offset) * RSD[4 + offset - 1];
      
      RSD[5 + offset - 1] /= BA(5,5,offset);
      RSD[4 + offset - 1] = (RSD[4 + offset - 1] - BA(4,5,offset) * RSD[5 + offset - 1]) / BA(4,4,offset);
      RSD[3 + offset - 1] = (RSD[3 + offset - 1] - BA(3,4,offset) * RSD[4 + offset - 1] - BA(3,5,offset) * RSD[5 + offset - 1]) / BA(3,3,offset);
      RSD[2 + offset - 1] = (RSD[2 + offset - 1] - BA(2,3,offset) * RSD[3 + offset - 1] - BA(2,4,offset) * RSD[4 + offset - 1] - BA(2,5,offset) * RSD[5 + offset - 1]) / BA(2,2,offset);
      RSD[1 + offset - 1] = (RSD[1 + offset - 1] - BA(1,2,offset) * RSD[2 + offset - 1] - BA(1,3,offset) * RSD[3 + offset - 1] - BA(1,4,offset) * RSD[4 + offset - 1] - BA(1,5,offset) * RSD[5 + offset - 1]) / BA(1,1,offset);
      
      for (j = ny - 1; j >= 1; --j) {
	YO(y, j)
        offset = 5 * (x + y + z);
	for (m = 1; m <= 5; ++m)
	  RSD[m + offset - 1] = RSD[m + offset - 1] - CA(m,1,offset) * RSD[1 + 5 * (x + Y_A[j] + z) - 1] - CA(m,2,offset) * RSD[2 + 5 * (x + Y_A[j] + z) - 1] - CA(m,3,offset) * RSD[3 + 5 * (x + Y_A[j] + z) - 1] - CA(m,4,offset) *
                    RSD[4 + 5 * (x + Y_A[j] + z) - 1] - 
                    CA(m,5,offset) * RSD[5 + 5 * (x + Y_A[j] + z) - 1];
	
	RSD[2 + offset - 1] -= BA(2,1,offset) * RSD[1 + offset - 1];
	RSD[3 + offset - 1] = RSD[3 + offset - 1] - BA(3,1,offset) * RSD[1 + offset - 1] - BA(3,2,offset) * RSD[2 + offset - 1];
	RSD[4 + offset - 1] = RSD[4 + offset - 1] - BA(4,1,offset) * RSD[1 + offset - 1] - BA(4,2,offset) * RSD[2 + offset - 1] - BA(4,3,offset) * RSD[3 + offset - 1];
	RSD[5 + offset - 1] = RSD[5 + offset - 1] - BA(5,1,offset) * RSD[1 + offset - 1] - BA(5,2,offset) * RSD[2 + offset - 1] - BA(5,3,offset) * RSD[3 + offset - 1] - BA(5,4,offset) * RSD[4 + offset - 1];
	
	RSD[5 + offset - 1] /= BA(5,5,offset);
                    RSD[4 + offset - 1] = (RSD[4 + offset - 1] - BA(4,5,offset) * RSD[5 + offset - 1]) / BA(4,4,offset);
	RSD[3 + offset - 1] = (RSD[3 + offset - 1] - BA(3,4,offset) * RSD[4 + offset - 1] - BA(3,5,offset) * RSD[5 + offset - 1]) / BA(3,3,offset);
	RSD[2 + offset - 1] = (RSD[2 + offset - 1] - BA(2,3,offset) * RSD[3 + offset - 1] - BA(2,4,offset) * RSD[4 + offset - 1] - BA(2,5,offset) * RSD[5 + offset - 1]) / BA(2,2,offset);
	RSD[1 + offset - 1] = (RSD[1 + offset - 1] - BA(1,2,offset) * RSD[2 + offset - 1] - BA(1,3,offset) * RSD[3 + offset - 1] - BA(1,4,offset) * RSD[4 + offset - 1] - BA(1,5,offset) * RSD[5 + offset - 1]) / BA(1,1,offset);
      }
    }
  }
  
  return(0);
} /* btridy_ */


int btridz_(nx, ny, nz)
int nx, ny, nz;
{
  int i, j, k, m, x, y, z, offset;
  double tm[25], tv[5], tmp;
  
/****solution of multiple, independent systems of block tridiagonal system
    s*/
/*   using Gaussian elimination (without pivoting) algorithm */
/*   (block size = 5 X 5) */
  
/* ***forward elimination phase */
  
  for (j = 2; j <= ny - 1; ++j) {
    YO(y, j)
    for (i = 2; i <= nx - 1; ++i) {
      XO(x, i)
      offset = 5 * (x + y + Z_A[0]);
      tmp = 1. / BA(1,1,offset);
      
      BA(2,1,offset) *= tmp;
      BA(3,1,offset) *= tmp;
      BA(4,1,offset) *= tmp;
      BA(5,1,offset) *= tmp;
      
      BA(2,2,offset) -= BA(2,1,offset) * BA(1,2,offset);
      
      tmp = 1. / BA(2,2,offset);
      
      BA(3,2,offset) = (BA(3,2,offset) - BA(3,1,offset) * BA(1,2,offset)) * tmp;
      BA(4,2,offset) = (BA(4,2,offset) - BA(4,1,offset) * BA(1,2,offset)) * tmp;
      BA(5,2,offset) = (BA(5,2,offset) - BA(5,1,offset) * BA(1,2,offset)) * tmp;
      
      BA(2,3,offset) -= BA(2,1,offset) * BA(1,3,offset);
      BA(3,3,offset) = BA(3,3,offset) - BA(3,1,offset) * BA(1,3,offset) - BA(3,2,offset) * BA(2,3,offset);
      
      tmp = 1. / BA(3,3,offset);
      
      BA(4,3,offset) = (BA(4,3,offset) - BA(4,1,offset) * BA(1,3,offset) - BA(4,2,offset) * BA(2,3,offset)) * tmp;
      BA(5,3,offset) = (BA(5,3,offset) - BA(5,1,offset) * BA(1,3,offset) - BA(5,2,offset) * BA(2,3,offset)) * tmp;
      
      BA(2,4,offset) -= BA(2,1,offset) * BA(1,4,offset);
      BA(3,4,offset) = BA(3,4,offset) - BA(3,1,offset) * BA(1,4,offset) - BA(3,2,offset) * BA(2,4,offset);
      BA(4,4,offset) = BA(4,4,offset) - BA(4,1,offset) * BA(1,4,offset) - BA(4,2,offset) * BA(2,4,offset) - BA(4,3,offset) * BA(3,4,offset);
      
      BA(5,4,offset) = (BA(5,4,offset) - BA(5,1,offset) * BA(1,4,offset) - BA(5,2,offset) * BA(2,4,offset) - BA(5,3,offset) * BA(3,4,offset)) /
                    BA(4,4,offset);
      
      BA(2,5,offset) -= BA(2,1,offset) * BA(1,5,offset);
      BA(3,5,offset) = BA(3,5,offset) - BA(3,1,offset) * BA(1,5,offset) - BA(3,2,offset) * BA(2,5,offset);
      BA(4,5,offset) = BA(4,5,offset) - BA(4,1,offset) * BA(1,5,offset) - BA(4,2,offset) * BA(2,5,offset) - BA(4,3,offset) * BA(3,5,offset);
                    
      BA(5,5,offset) = BA(5,5,offset) - BA(5,1,offset) * BA(1,5,offset) - BA(5,2,offset) * BA(2,5,offset) - BA(5,3,offset) * BA(3,5,offset) - 
                    BA(5,4,offset) * BA(4,5,offset);
      
      for (k = 2; k <= nz; ++k) {
	ZO(z, k-1)
	offset = 5 * (x + y + z);
	for (m = 1; m <= 5; ++m) {
	  
	  tm[m * 5 - 5] = CA(1,m,offset);
	  tm[m * 5 - 4] = CA(2,m,offset) - BA(2,1,offset) * tm[m * 5 - 5];
	  tm[m * 5 - 3] = CA(3,m,offset) - BA(3,1,offset) * tm[m * 5 - 
                    5] - BA(3,2,offset) * tm[m * 5 - 4];
	  tm[m * 5 - 2] = CA(4,m,offset) - BA(4,1,offset) * tm[m * 5 - 
                    5] - BA(4,2,offset) * tm[m * 5 - 4] - BA(4,3,offset) * tm[m * 
                    5 - 3];
	  tm[m * 5 - 1] = CA(5,m,offset) - BA(5,1,offset) * tm[m * 5 - 
                    5] - BA(5,2,offset) * tm[m * 5 - 4] - BA(5,3,offset) * tm[m * 
                    5 - 3] - BA(5,4,offset) * tm[m * 5 - 2];
	}
	
	for (m = 1; m <= 5; ++m) {
	  
	  tm[m * 5 - 1] /= BA(5,5,offset);
	  tm[m * 5 - 2] = (tm[m * 5 - 2] - BA(4,5,offset) * tm[m * 5 - 
                    1]) / BA(4,4,offset);
	  tm[m * 5 - 3] = (tm[m * 5 - 3] - BA(3,4,offset) * tm[m * 5 - 
                    2] - BA(3,5,offset) * tm[m * 5 - 1]) / BA(3,3,offset);
	  tm[m * 5 - 4] = (tm[m * 5 - 4] - BA(2,3,offset) * tm[m * 5 - 
                    3] - BA(2,4,offset) * tm[m * 5 - 2] - BA(2,5,offset) * tm[m * 
                    5 - 1]) / BA(2,2,offset);
	  tm[m * 5 - 5] = (tm[m * 5 - 5] - BA(1,2,offset) * tm[m * 5 - 
                    4] - BA(1,3,offset) * tm[m * 5 - 3] - BA(1,4,offset) * tm[m * 
                    5 - 2] - BA(1,5,offset) * tm[m * 5 - 1]) / BA(1,1,offset);
	}
	
	ZO(z, k)
	offset = 5 * (x + y + z);
	for (m = 1; m <= 5; ++m) {

	  BA(1,m,offset) = BA(1,m,offset) 
                    - AA(1,1,offset) * tm[m * 5 - 5] - AA(1,2,offset) * tm[m * 5 - 4] - AA(1,3,offset) * 
                    tm[m * 5 - 3] - AA(1,4,offset) * tm[m * 5 - 2] - AA(1,5,offset) * tm[m * 5 - 1];
	  
	  BA(2,m,offset) = BA(2,m,offset) 
                    - AA(2,1,offset) * tm[m * 5 - 5] - AA(2,2,offset) * tm[m * 5 - 4] - AA(2,3,offset) * 
                    tm[m * 5 - 3] - AA(2,4,offset) * tm[m * 5 - 2] - AA(2,5,offset) * tm[m * 5 - 1];
	  
	  BA(3,m,offset) = BA(3,m,offset) 
                    - AA(3,1,offset) * tm[m * 5 - 5] - AA(3,2,offset) * tm[m * 5 - 4] - AA(3,3,offset) * 
                    tm[m * 5 - 3] - AA(3,4,offset) * tm[m * 5 - 2] - AA(3,5,offset) * tm[m * 5 - 1];
	  
	  BA(4,m,offset) = BA(4,m,offset) 
                    - AA(4,1,offset) * tm[m * 5 - 5] - AA(4,2,offset) * tm[m * 5 - 4] - AA(4,3,offset) * 
                    tm[m * 5 - 3] - AA(4,4,offset) * tm[m * 5 - 2] - AA(4,5,offset) * tm[m * 5 - 1];
	  
	  BA(5,m,offset) = BA(5,m,offset) 
                    - AA(5,1,offset) * tm[m * 5 - 5] - AA(5,2,offset) * tm[m * 5 - 4] - AA(5,3,offset) * 
                    tm[m * 5 - 3] - AA(5,4,offset) * tm[m * 5 - 2] - AA(5,5,offset) * tm[m * 5 - 1];
	}
	offset = 5 * (x + y + Z_A[k-2]);
	tv[0] = RSD[1 + offset - 1];
	tv[1] = RSD[2 + offset - 1] - BA(2,1,offset) * tv[0];
	tv[2] = RSD[3 + offset - 1] - BA(3,1,offset) * tv[0] - BA(3,2,offset) * tv[1];
	tv[3] = RSD[4 + offset - 1] - BA(4,1,offset) * tv[0] - BA(4,2,offset) * tv[1] - BA(4,3,offset) * tv[2];
	tv[4] = RSD[5 + offset - 1] - BA(5,1,offset) * tv[0] - BA(5,2,offset) * tv[1] - BA(5,3,offset) * tv[2] - BA(5,4,offset) * 
                    tv[3];
	
	tv[4] /= BA(5,5,offset);
	tv[3] = (tv[3] - BA(4,5,offset) * tv[4]) / BA(4,4,offset);
	tv[2] = (tv[2] - BA(3,4,offset) * tv[3] - BA(3,5,offset) * tv[4]) / BA(3,3,offset);
	tv[1] = (tv[1] - BA(2,3,offset) * tv[2] - BA(2,4,offset) * tv[3] - BA(2,5,offset) * tv[4]) 
                    / BA(2,2,offset);
	tv[0] = (tv[0] - BA(1,2,offset) * tv[1] - BA(1,3,offset) * tv[2] - BA(1,4,offset) * tv[3] 
                    - BA(1,5,offset) * tv[4]) / BA(1,1,offset);
	
        ZO(z, k)
	offset = 5 * (x + y + z);
	for (m = 1; m <= 5; ++m)
	  RSD[m + offset - 1] = RSD[m + offset - 1] - AA(m,1,offset) * tv[0] - AA(m,2,offset) *
                    tv[1] - AA(m,3,offset) * tv[2] - AA(m,4,offset) * tv[3] - AA(m,5,offset) * tv[4];
	
	tmp = 1. / BA(1,1,offset);
	
	BA(2,1,offset) *= tmp;
	BA(3,1,offset) *= tmp;
	BA(4,1,offset) *= tmp;
	BA(5,1,offset) *= tmp;
	
	BA(2,2,offset) -= BA(2,1,offset) * BA(1,2,offset);
	
	tmp = 1. / BA(2,2,offset);
	
	BA(3,2,offset) = (BA(3,2,offset) - BA(3,1,offset) * BA(1,2,offset)) * tmp;
	BA(4,2,offset) = (BA(4,2,offset) - BA(4,1,offset) * BA(1,2,offset)) * tmp;
	BA(5,2,offset) = (BA(5,2,offset) - BA(5,1,offset) * BA(1,2,offset)) * tmp;
	
	BA(2,3,offset) -= BA(2,1,offset) * BA(1,3,offset);
	BA(3,3,offset) = BA(3,3,offset) - BA(3,1,offset) * BA(1,3,offset) - BA(3,2,offset) * 
                    BA(2,3,offset);
	
	tmp = 1. / BA(3,3,offset);
	
	BA(4,3,offset) = (BA(4,3,offset) - BA(4,1,offset) * BA(1,3,offset) - BA(4,2,offset) * 
                    BA(2,3,offset)) 
                    * tmp;
	BA(5,3,offset) = (BA(5,3,offset) - BA(5,1,offset) * BA(1,3,offset) - BA(5,2,offset) * 
                    BA(2,3,offset)) 
                    * tmp;
	
	BA(2,4,offset) -= BA(2,1,offset) * BA(1,4,offset);
	BA(3,4,offset) = BA(3,4,offset) - BA(3,1,offset) * BA(1,4,offset) - BA(3,2,offset) * 
                    BA(2,4,offset);
	BA(4,4,offset) = BA(4,4,offset) - BA(4,1,offset) * BA(1,4,offset) - BA(4,2,offset) * 
                    BA(2,4,offset) 
                    - BA(4,3,offset)
                    * BA(3,4,offset);
	
	BA(5,4,offset) = (BA(5,4,offset) - BA(5,1,offset) * BA(1,4,offset) - BA(5,2,offset) * 
                    BA(2,4,offset) 
                    - BA(5,3,offset)
                    * BA(3,4,offset)) / BA(4,4,offset);
	
	BA(2,5,offset) -= BA(2,1,offset) * BA(1,5,offset);
	BA(3,5,offset) = BA(3,5,offset) - BA(3,1,offset) * BA(1,5,offset) - BA(3,2,offset) * 
                    BA(2,5,offset);
	BA(4,5,offset) = BA(4,5,offset) - BA(4,1,offset) * BA(1,5,offset) - BA(4,2,offset) * 
                    BA(2,5,offset) 
                    - BA(4,3,offset)
                    * BA(3,5,offset);
	
	BA(5,5,offset) = BA(5,5,offset) - BA(5,1,offset) * BA(1,5,offset) - BA(5,2,offset) * 
                    BA(2,5,offset) 
                    - BA(5,3,offset)
                    * BA(3,5,offset) - BA(5,4,offset) * BA(4,5,offset);
      }
    }
  }
  
/* ***back-substitution phase */
  
  for (j = 2; j <= ny - 1; ++j) {
    YO(y, j)
    for (i = 2; i <= nx - 1; ++i) {
      XO(x, i)
      offset = 5 * (x + y + Z_A[nz-1]);
      RSD[2 + offset - 1] -= BA(2,1,offset) * RSD[1 + offset - 1];
      RSD[3 + offset - 1] = RSD[3 + offset - 1] - BA(3,1,offset) * RSD[1 + offset - 1] - BA(3,2,offset) * RSD[2 + offset - 1];
      RSD[4 + offset - 1] = RSD[4 + offset - 1] - BA(4,1,offset) * RSD[1 + offset - 1] - BA(4,2,offset) * RSD[2 + offset - 1] - BA(4,3,offset) * RSD[3 + offset - 1];
      RSD[5 + offset - 1] = RSD[5 + offset - 1] - BA(5,1,offset) * RSD[1 + offset - 1] - BA(5,2,offset) * RSD[2 + offset - 1] - BA(5,3,offset) * RSD[3 + offset - 1] - BA(5,4,offset) * RSD[4 + offset - 1];
      
      RSD[5 + offset - 1] /= BA(5,5,offset);
      RSD[4 + offset - 1] = (RSD[4 + offset - 1] - BA(4,5,offset) * RSD[5 + offset - 1]) / BA(4,4,offset);
      RSD[3 + offset - 1] = (RSD[3 + offset - 1] - BA(3,4,offset) * RSD[4 + offset - 1] - BA(3,5,offset) * RSD[5 + offset - 1]) / BA(3,3,offset);
      RSD[2 + offset - 1] = (RSD[2 + offset - 1] - BA(2,3,offset) * RSD[3 + offset - 1] - BA(2,4,offset) * RSD[4 + offset - 1] - BA(2,5,offset) * RSD[5 + offset - 1]) / BA(2,2,offset);
      RSD[1 + offset - 1] = (RSD[1 + offset - 1] - BA(1,2,offset) * RSD[2 + offset - 1] - BA(1,3,offset) * RSD[3 + offset - 1] - BA(1,4,offset) * RSD[4 + offset - 1] - BA(1,5,offset) * RSD[5 + offset - 1]) / BA(1,1,offset);
      
      for (k = nz - 1; k >= 1; --k) {
	ZO(z, k)
	offset = 5 * (x + y + z);
	for (m = 1; m <= 5; ++m)
	  RSD[m + offset - 1] = RSD[m + offset - 1] - CA(m,1,offset) * RSD[1 + 5 * (x + y + Z_A[k]) - 1] - CA(m,2,offset) * RSD[2 + 5 * (x + y + Z_A[k]) - 1] - CA(m,3,offset) *
                    RSD[3 + 5 * (x + y + Z_A[k]) - 1] 
                    - CA(m,4,offset) * RSD[4 + 5 * (x + y + Z_A[k]) - 1] - 
		      CA(m,5,offset) * RSD[5 + 5 * (x + y + Z_A[k]) - 1];
	
	RSD[2 + offset - 1] -= BA(2,1,offset) * RSD[1 + offset - 1];
	RSD[3 + offset - 1] = RSD[3 + offset - 1] - BA(3,1,offset) * RSD[1 + offset - 1] - BA(3,2,offset) * RSD[2 + offset - 1];
	RSD[4 + offset - 1] = RSD[4 + offset - 1] - BA(4,1,offset) * RSD[1 + offset - 1] - BA(4,2,offset) * RSD[2 + offset - 1] - BA(4,3,offset) * RSD[3 + offset - 1];
	RSD[5 + offset - 1] = RSD[5 + offset - 1] - BA(5,1,offset) * RSD[1 + offset - 1] - BA(5,2,offset) * RSD[2 + offset - 1] - BA(5,3,offset) * RSD[3 + offset - 1] - BA(5,4,offset) * RSD[4 + offset - 1] ;
	
	RSD[5 + offset - 1] /= BA(5,5,offset);
	RSD[4 + offset - 1] = (RSD[4 + offset - 1] - BA(4,5,offset) * RSD[5 + offset - 1]) / BA(4,4,offset);
	RSD[3 + offset - 1] = (RSD[3 + offset - 1] - BA(3,4,offset) * RSD[4 + offset - 1] - BA(3,5,offset) * RSD[5 + offset - 1]) / BA(3,3,offset);
	RSD[2 + offset - 1] = (RSD[2 + offset - 1] - BA(2,3,offset) * RSD[3 + offset - 1] - BA(2,4,offset) * RSD[4 + offset - 1] - BA(2,5,offset) * RSD[5 + offset - 1]) / BA(2,2,offset);
	RSD[1 + offset - 1] = (RSD[1 + offset - 1] - BA(1,2,offset) * RSD[2 + offset - 1] - BA(1,3,offset) * RSD[3 + offset - 1] - BA(1,4,offset) * RSD[4 + offset - 1] - BA(1,5,offset) * RSD[5 + offset - 1]) / BA(1,1,offset);
      }
    }
  }
  return(0);
} /* btridz_ */

