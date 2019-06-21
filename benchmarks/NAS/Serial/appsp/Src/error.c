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

/* error.c */

#define MAIN extern
#include "common.h"

#define imax(dim1) ONE_D_ARRAY(IMAX, dim1)
#define jmax(dim1) ONE_D_ARRAY(JMAX, dim1)
#define kmax(dim1) ONE_D_ARRAY(KMAX, dim1)
#define u000ijk(dim1) ONE_D_ARRAY(U000IJK, dim1)
#define errmax(dim1) ONE_D_ARRAY(ERRMAX, dim1)

/* compute the solution error */

error()
{
  int lnorm;
  int i,j,k,m;
  static double tmp;
  int IMAX[5], JMAX[5], KMAX[5];
  static double U000IJK[5], ERRMAX[5];
  
  lnorm = 2;
  if ( lnorm == 1 )
    {
      for (m=1;m<=5;m++) 
	{
	  errmax(m) = - 1.0e+20;
	}

      for (k=2;k<=nz-1;k++)
	{
	  for (j=2;j<=ny-1;j++)
	    {
	      for (i=2;i<=nx-1;i++)
		{
		  exact ( i, j, k, &u000ijk(1) );
		  for (m=1;m<=5;m++) 
		    {
		      tmp = fabs ( u000ijk(m) - u(m,i,j,k) );
		      if ( tmp == errmax(m) ) 
			{
			  errmax(m) = tmp;
			  imax(m) = i;
			  jmax(m) = j;
			  kmax(m) = k;
			}
		    }
		}
            }
	}

      fprintf(outfp, "     max. error in soln. to first pde  = %e\n");
      fprintf(outfp, "     and its location                  = (%4d,%4d,%4d)\n",  
	      errmax(1), imax(1), jmax(1), kmax(1));

      fprintf(outfp, "     max. error in soln. to second pde = %e\n");
      fprintf(outfp, "     and its location                  = (%4d,%4d,%4d)\n",  
	      errmax(2), imax(2), jmax(2), kmax(2));

      fprintf(outfp, "     max. error in soln. to third pde  = %e\n");
      fprintf(outfp, "     and its location                  = (%4d,%4d,%4d)\n",  
	      errmax(3), imax(3), jmax(3), kmax(3));
      fprintf(outfp, "     max. error in soln. to fourth pde = %e\n");
      fprintf(outfp, "     and its location                  = (%4d,%4d,%4d)\n",  
	      errmax(4), imax(4), jmax(4), kmax(4));
      fprintf(outfp, "     max. error in soln. to fifth pde  = %e\n");
      fprintf(outfp, "     and its location                  = (%4d,%4d,%4d)\n",  
	      errmax(5), imax(5), jmax(5), kmax(5));
      fflush(outfp);
    }
  else if ( lnorm == 2 ) 
    {
      for (m=1;m<=5;m++)
	{
	  errnm(m) = 0.0;
	}

      for (k=2;k<=nz-1;k++) 
	{
	  for (j=2;j<=ny-1;j++) 
	    {
	      for (i=2;i<=nx-1;i++) 
		{
		  exact ( i, j, k, &u000ijk(1) );
		  for (m=1;m<=5;m++) 
		    {
		      tmp = ( u000ijk(m) - u(m,i,j,k) );
		      errnm(m) = errnm(m) + tmp*tmp;
		    }
		}
            }
	}
      for (m=1;m<=5;m++) 
	{
	  errnm(m) = sqrt ( errnm(m) / ( (nx-2)*(ny-2)*(nz-2) ) );
	}

      fprintf(outfp, " RMS-norm of error in soln. to first pde  = %e\n", errnm(1));
      fprintf(outfp, " RMS-norm of error in soln. to second pde = %e\n", errnm(2));
      fprintf(outfp, " RMS-norm of error in soln. to third pde  = %e\n", errnm(3));
      fprintf(outfp, " RMS-norm of error in soln. to fourth pde = %e\n", errnm(4));
      fprintf(outfp, " RMS-norm of error in soln. to fifth pde  = %e\n\n", errnm(5));
      fflush(outfp);
    }
}

