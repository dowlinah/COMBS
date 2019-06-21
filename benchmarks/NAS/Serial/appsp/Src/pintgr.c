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

/* pintgr.c */

#define MAIN extern
#include "common.h"

#define PHI_SIZE		ISIZ1*ISIZ2
#define phi1(dim1, dim2) 	(*(PHI1 + (dim2 - 1) * ISIZ2 + (dim1 -1 )))
#define phi2(dim1, dim2) 	(*(PHI2 + (dim2 - 1) * ISIZ2 + (dim1 -1 )))

/* compute the surface integral */

pintgr()
{
  double PHI1[PHI_SIZE], PHI2[PHI_SIZE];
  int i,j,k;
  double frc1, frc2, frc3;
  
  for (j=ji1;j<=ji2;j++)
    {
      for (i=ii1;i<=ii2;i++)
	{
	  phi1(i,j) = c2*(  u(5,i,j,ki1)
			  - 0.50 * (  u(2,i,j,ki1)*u(2,i,j,ki1)
				    + u(3,i,j,ki1)*u(3,i,j,ki1)
				    + u(4,i,j,ki1)*u(4,i,j,ki1))
			  / u(1,i,j,ki1) );
	  phi2(i,j) = c2*(u(5,i,j,ki2)
			  - 0.50 * (u(2,i,j,ki2)*u(2,i,j,ki2)
				    + u(3,i,j,ki2)*u(3,i,j,ki2)
				    + u(4,i,j,ki2)*u(4,i,j,ki2) )
			  / u(1,i,j,ki2) );
	}
    }
  
  frc1 = 0.0;
  for (j=ji1;j<=ji2-1;j++)
    {
      for (i=ii1;i<=ii2-1;i++)
	{
	  frc1 = frc1 + (  phi1(i,j)
			 + phi1(i+1,j)
			 + phi1(i,j+1)
			 + phi1(i+1,j+1)
			 + phi2(i,j)
			 + phi2(i+1,j)
			 + phi2(i,j+1)
			 + phi2(i+1,j+1) );
	}
    }
  frc1 = dxi * deta * frc1;
  for (k=ki1;k<=ki2;k++)
    {
      for (i=ii1;i<=ii2;i++)
	{
	  phi1(i,k) = c2*(  u(5,i,ji1,k)
			  - 0.50 * ( u(2,i,ji1,k)*u(2,i,ji1,k)
				    + u(3,i,ji1,k)*u(3,i,ji1,k)
				    + u(4,i,ji1,k)*u(4,i,ji1,k) )
			  / u(1,i,ji1,k) );
	  phi2(i,k) = c2*(  u(5,i,ji2,k)
			  - 0.50 * (u(2,i,ji2,k)*u(2,i,ji2,k)
				    + u(3,i,ji2,k)*u(3,i,ji2,k)
				    + u(4,i,ji2,k)*u(4,i,ji2,k) )
			  / u(1,i,ji2,k) );
	}
    }
  frc2 = 0.0;
  for(k=ki1;k<=ki2-1;k++)
    {
      for (i=ii1;i<=ii2-1;i++)
	{
	  frc2 = frc2 + (  phi1(i,k)
			 + phi1(i+1,k)
			 + phi1(i,k+1)
			 + phi1(i+1,k+1)
			 + phi2(i,k)
			 + phi2(i+1,k)
			 + phi2(i,k+1)
			 + phi2(i+1,k+1) );
	}
    }
  frc2 = dxi * dzeta * frc2;
  for (k=ki1;k<=ki2;k++)
    {
      for (j=ji1;j<=ji2;j++)
	{
	  phi1(j,k) = c2*(  u(5,ii1,j,k)
			  - 0.50 * (u(2,ii1,j,k)*u(2,ii1,j,k)
				    + u(3,ii1,j,k)*u(3,ii1,j,k)
				    + u(4,ii1,j,k)*u(4,ii1,j,k) )
			  / u(1,ii1,j,k) );
	  phi2(j,k) = c2*(  u(5,ii2,j,k)
			  - 0.50 * (  u(2,ii2,j,k)*u(2,ii2,j,k)
				    + u(3,ii2,j,k)*u(3,ii2,j,k)
				    + u(4,ii2,j,k)*u(4,ii2,j,k) )
			  / u(1,ii2,j,k) );
	}
    }
  frc3 = 0.0;
  for (k=ki1;k<=ki2-1;k++)
    {
      for (j=ji1;j<=ji2-1;j++)
	{
	  frc3 = frc3 + (  phi1(j,k)
			 + phi1(j+1,k)
			 + phi1(j,k+1)
			 + phi1(j+1,k+1)
			 + phi2(j,k)
			 + phi2(j+1,k)
			 + phi2(j,k+1)
			 + phi2(j+1,k+1) );
	}
    }
  frc3 = deta * dzeta * frc3;
  frc = 0.25 * ( frc1 + frc2 + frc3 );
  fprintf(outfp, "\n     surface integral = %e\n\n\n\n", frc);
}
