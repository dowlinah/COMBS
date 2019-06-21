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

/* setbv.c */

#define MAIN extern
#include "common.h"

/* set the boundary values of dependent values */

setbv()
{
  int i, j, k;

  /* set the dependent variable values along the top and bottom faces */

  for (j=1;j<=ny;j++) 
    {
      for (i=1;i<=nx;i++) 
	{
	  exact ( i, j, 1, &u( 1, i, j, 1 ) );
	  exact ( i, j, nz,&u( 1, i, j, nz ) );
	}
    }

  /* set the dependent variable values along north and south faces */

  for (k=1;k<=nz;k++) 
    {
      for (i=1;i<=nx;i++) 
	{
      	  exact ( i, 1, k, &u( 1, i, 1, k ) );
	  exact ( i, ny, k,&u( 1, i, ny, k ) );
	}
    }

  /* set the dependent variable values along east and west faces */

  for (k=1;k<=nz;k++) 
    {
      for (j=1;j<=ny;j++)
	{
	  exact ( 1, j, k, &u( 1, 1, j, k ) );
	  exact ( nx, j, k, &u( 1, nx, j, k ) );
	}
    }
}
      
      
