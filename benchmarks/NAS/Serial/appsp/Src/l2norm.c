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

/* l2norm.c */

#define MAIN extern
#include "common.h"

#define v(dim1, dim2, dim3, dim4)	FOUR_D_ARRAY(V, dim1, dim2, dim3, dim4)
#define sum(dim1)			ONE_D_ARRAY(SUM, dim1)

/* compute the l2-norm of vector v */

l2norm ( ldx, ldy, ldz, nx, ny, nz, V, SUM )
int ldx, ldy, ldz;
int nx, ny, nz;
double *V;
double *SUM;
{
  int i,j,k,m;

  for (m=1;m<=5;m++)
    {
      sum(m) = 0.0;
    }
  for (k=2;k<=nz-1;k++) 
    {
      for (j=2;j<=ny-1;j++)
	{
	  for (i=2;i<=nx-1;i++)
	    {
	      for (m=1;m<=5;m++) 
		{
                  sum(m) = sum(m) + v(m,i,j,k) * v(m,i,j,k);
		}
	    }
	}
    }

  for (m=1;m<=5;m++)
    {
      sum(m) = sqrt ( sum(m) / ( (nx-2) * (ny-2) * (nz-2) ) );
    }
}

