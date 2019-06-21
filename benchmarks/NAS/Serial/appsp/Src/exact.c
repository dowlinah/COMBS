/*
|| NAS CFD code: APPSP
||
|| Fortran version:
||
|| Author: Sisira Weeratunga
||         NASA Ames Research Center
||         (10/25/90)
||
|| This is the C version of the Fortran code developed at NASA Ames Research
|| Center.  We converted the code to C as a part of our class project for 
|| CS 838-3/ChE 562  offered in Spring 1993 by Mark D. Hill, Sangtae Kim and
|| Mary Vernon at the University of Wisconsin at Madison.  CS 838-3/ChE 562 
|| was an experimental course that brought computer scientists and computation
|| scientists together to promote interdisciplinary research.
||
||					June 18, 1993.
||
||					Shubhendu S. Mukherjee (Shubu)
||					shubu@cs.wisc.edu
||					Computer Sciences Department
||
||					Iasonas Moustakis
||					iasonas@luther.che.wisc.edu
||					Chemical Engineering Department
||					
||					University of Wisconsin at Madison
||
*/

/* exact.c */

#define MAIN extern
#include "common.h"

#define u000ijk(dim1) 	ONE_D_ARRAY(U000IJK, dim1)

/* compute the exact solution at (i,j,k) */

exact(i,j,k,U000IJK)
int i,j,k;
double *U000IJK;
{
  double xi, eta, zeta;
  int m;

  /* compute the exact solution at (i,j,k) */
  
  xi  = (( double )  ( i - 1 ))  / ( nx - 1 );
  eta  = (( double ) ( j - 1 ))  / ( ny - 1 );
  zeta = (( double ) ( k - 1 ))  / ( nz - 1 );
  
  for(m = 1;m<= 5;m++)
    {
      u000ijk(m) =  ce(m,1) 
	+ ce(m,2) * xi 
	+ ce(m,3) * eta
	+ ce(m,4) * zeta
	+ ce(m,5) * xi * xi
	+ ce(m,6) * eta * eta
        + ce(m,7) * zeta * zeta
	+ ce(m,8) * xi * xi * xi
	+ ce(m,9) * eta * eta * eta
	+ ce(m,10) * zeta * zeta * zeta
	+ ce(m,11) * xi * xi * xi * xi
	+ ce(m,12) * eta * eta * eta * eta
	+ ce(m,13) * zeta * zeta * zeta * zeta;
    }
}
