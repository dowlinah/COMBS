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

/* setiv.c */

#define MAIN extern
#include "common.h"

/*
|| set the initial values of independent variables based on tri-linear
|| interpolation of boundary values in the computational space
*/

setiv()
{
  int i, j, k, m;
  double xi, eta, zeta, pxi, peta, pzeta;
  
  for (k = 2; k <= nz-1; k ++)
    {
      zeta = (( double) (k-1)) / (nz-1);
      for (j = 2; j <= ny-1; j ++)
	{
	  eta = (( double )(j-1)) / (ny-1);
	  for (i = 2; i <= nx-1; i ++)
	    {
	      xi = (( double ) (i-1))  / (nx-1);
	      
	      for (m=1;m<=5;m++)
		{
		  pxi =   ( 1.0 - xi ) * u(m,1,j,k)
		    + xi   * u(m,nx,j,k);
		  
		  peta =  ( 1.0 - eta ) * u(m,i,1,k)
		    + eta   * u(m,i,ny,k);
		  
		  pzeta = ( 1.0 - zeta ) * u(m,i,j,1)
		    + zeta   * u(m,i,j,nz);
		  
		  u( m, i, j, k ) = pxi + peta + pzeta
		    - pxi * peta - peta * pzeta - pzeta * pxi
		      + pxi * peta * pzeta;
		}
	    }
	}
    }
}
