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

/* tinv.c */

#define MAIN extern
#include "common.h"

static double ru1, uu, vv, ww, q, ac, ac2, r1, r2,
              r3, r4, r5, t1, t2, t3, alph;

/* block-diagonal matrix-vector multiplication */

txinvr()
{
  double bt = sqrt ( 0.50 );
  int i,j,k;
  
  for (k=2;k<=nz-1;k++)
    {
      for (j=2;j<=ny-1;j++)
	{
	  for (i=2;i<=nx-1;i++)
	    {
	      ru1 = 1.0 / u(1,i,j,k);	
	      uu = ru1 * u(2,i,j,k);
	      vv = ru1 * u(3,i,j,k);
	      ww = ru1 * u(4,i,j,k);
	      q = 0.50 * (   uu*uu
			  +  vv*vv
			  +  ww*ww );
	      ac2 = c1 * c2 * ( ru1 * u(5,i,j,k) - q );
	      if ( ac2 <= 0.0 )
		{
		  fprintf(outfp, "Speed of sound is zero");
		  exit();
		}
	      ac = sqrt ( ac2 );
	      r1 = rsd(1,i,j,k);
	      r2 = rsd(2,i,j,k);
	      r3 = rsd(3,i,j,k);
	      r4 = rsd(4,i,j,k);
	      r5 = rsd(5,i,j,k);
	      t1 = ( c2 / ac2 ) * ( q * r1 - uu * r2
				   - vv * r3
				   - ww * r4
				   +      r5 );
	      rsd(1,i,j,k) = r1 - t1 ;
	      rsd(2,i,j,k) = - ru1 * ( ww * r1 - r4 );
	      rsd(3,i,j,k) =   ru1 * ( vv * r1 - r3 );
	      t2 = bt * ru1 * ( uu * r1 - r2 );
	      t3 = ( bt * ru1 * ac ) * t1;
	      rsd(4,i,j,k) = - t2 + t3;
	      rsd(5,i,j,k) =   t2 + t3;
	    }
	}
    }
}

/* block-diagonal matrix-vector multipication */

tzetar()
{
  double bt = sqrt ( 0.50 );
  int i, j, k;
  
  for (k=2;k<=nz-1;k++)
    {
      for (j=2;j<=ny-1;j++)
	{
	  for (i=2;i<=nx-1;i++)
	    {
	      ru1 = 1.0 / u(1,i,j,k);
	      uu = ru1 * u(2,i,j,k);
	      vv = ru1 * u(3,i,j,k);
	      ww = ru1 * u(4,i,j,k);
	      q = 0.50 * (   uu*uu
			  +  vv*vv
			  +  ww*ww );
	      ac2 = c1 * c2 * ( ru1 * u(5,i,j,k) - q );
	      ac = sqrt ( ac2 );
	      alph = ( bt * u(1,i,j,k) ) / ac;
	      r1 = rsd(1,i,j,k);
	      r2 = rsd(2,i,j,k);
	      r3 = rsd(3,i,j,k);
	      r4 = rsd(4,i,j,k);
	      r5 = rsd(5,i,j,k);
	      t1 = alph * ( r4 + r5 );
	      t2 = r3 + t1;
	      t3 = bt * u(1,i,j,k) * ( r4 - r5 );
	      rsd(1,i,j,k) = t2;
	      rsd(2,i,j,k) = - u(1,i,j,k) * r2 + uu * t2;
	      rsd(3,i,j,k) =   u(1,i,j,k) * r1 + vv * t2;
	      rsd(4,i,j,k) =   ww * t2 + t3;
	      rsd(5,i,j,k) =  u(1,i,j,k) * ( - uu * r2 + vv * r1 )
		+ q * t2
		  + ( ac2 / c2 ) * t1
		    + ww * t3;
	    }
	}
    }
}


