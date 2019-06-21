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

/* invr.c */

#define MAIN extern
#include "common.h"

#define PHI_SIZE		ISIZ1*ISIZ2
#define phi1(dim1, dim2) 	(*(PHI1 + (dim2 - 1) * ISIZ2 + (dim1 -1 )))
#define phi2(dim1, dim2) 	(*(PHI2 + (dim2 - 1) * ISIZ2 + (dim1 -1 )))

/* block-diagonal matrix-vector multiply */

ninvr()
{
  double bt = sqrt(0.5);
  int i,j,k;
  double r1, r2, r3, r4, r5, t1, t2;
  
  
  for (k=2;k<=nz-1;k++)
    {
      for (j=2;j<=ny-1;j++)
	{
	  for (i=2;i<=nx-1;i++)
	    {
	      r1 = rsd(1,i,j,k);
	      r2 = rsd(2,i,j,k);
	      r3 = rsd(3,i,j,k);
	      r4 = rsd(4,i,j,k);
	      r5 = rsd(5,i,j,k);
	      rsd(1,i,j,k) = - r2;
	      rsd(2,i,j,k) =   r1;
	      rsd(3,i,j,k) = bt * ( r4 - r5 );
	      t1 = bt * r3;
	      t2 = 0.50 * ( r4 + r5 );
	      rsd(4,i,j,k) = - t1 + t2;
	      rsd(5,i,j,k) =   t1 + t2;
	    }
	}
    }
}

/* block-diagonal matrix-vector multiply */

pinvr()
{
  double  bt = sqrt ( 0.50 );
  int i,j,k;
  double r1, r2, r3, r4, r5, t1, t2;  
  
  for (k=2;k<=nz-1;k++)
    {
      for (j=2;j<=ny-1;j++)
	{
	  for (i=2;i<=nx-1;i++)
	    {
	      r1 = rsd(1,i,j,k);
	      r2 = rsd(2,i,j,k);
	      r3 = rsd(3,i,j,k);
	      r4 = rsd(4,i,j,k);
	      r5 = rsd(5,i,j,k);
	      rsd(1,i,j,k) = bt * ( r4 - r5 );
	      rsd(2,i,j,k) = - r3;
	      rsd(3,i,j,k) = r2;
	      t1 = bt * r1;
	      t2 = 0.50 * ( r4 + r5 );
	      rsd(4,i,j,k) = - t1 + t2;
	      rsd(5,i,j,k) =   t1 + t2;
	    }
	}
    }
}	





