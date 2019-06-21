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

/* spenta.c */

#define MAIN extern
#include "common.h"

#define f(dim1, dim2, dim3, dim4)	FOUR_D_ARRAY(F, dim1, dim2, dim3, dim4)

/* temporary variables */
int i,j,k;
static double tmp1, tmp2;

/* solution of multiple, independent systems of penta-diagonal systems */

spentax ( isiz1, isiz2, isiz3, m, nx, ny, nz, A, B, C, D, E, F )
int isiz1, isiz2, isiz3;
int m, nx,ny, nz;
double *A, *B, *C, *D, *E, *F;
{
  /* forward elimination */
  
  for (k=2;k<=nz-1;k++)
    {
      for (j=2;j<=ny-1;j++)
	{
	  for (i=2;i<=nx-1;i++)
	    {
	      tmp1 = b(i,j,k) / c(i-1,j,k);
	      tmp2 = a(i+1,j,k) / c(i-1,j,k);
	      c(i,j,k) = c(i,j,k) - tmp1 * d(i-1,j,k);
	      d(i,j,k) = d(i,j,k) - tmp1 * e(i-1,j,k) ;
	      f(m,i,j,k) = f(m,i,j,k) - tmp1 * f(m,i-1,j,k);
	      b(i+1,j,k) = b(i+1,j,k) - tmp2 * d(i-1,j,k);
	      c(i+1,j,k) = c(i+1,j,k) - tmp2 * e(i-1,j,k);
	      f(m,i+1,j,k) = f(m,i+1,j,k) - tmp2 * f(m,i-1,j,k);
	    }
	  tmp1 = b(nx,j,k) / c(nx-1,j,k);
	  c(nx,j,k) = c(nx,j,k) - tmp1 * d(nx-1,j,k);
	  d(nx,j,k) = d(nx,j,k) - tmp1 * e(nx-1,j,k);
	  f(m,nx,j,k) = f(m,nx,j,k) - tmp1 * f(m,nx-1,j,k);
	}
    }

  /* back-substitution phase */

  for (k=2;k<=nz-1;k++)
    {
      for (j=2;j<=ny-1;j++)
	{
	  f(m,nx,j,k) = f(m,nx,j,k) / c(nx,j,k);
	  f(m,nx-1,j,k) = ( f(m,nx-1,j,k) - d(nx-1,j,k)*f(m,nx,j,k) )
	    / c(nx-1,j,k);
	  for (i=nx-2;i>=1;i--)
	    {
	      f(m,i,j,k) = ( f(m,i,j,k) - d(i,j,k)*f(m,i+1,j,k)
			    - e(i,j,k)*f(m,i+2,j,k) ) / c(i,j,k);
	    }
	}
    }
}

/* solution of multiple, independent systems of penta-diagonal systems */

spentax3 ( isiz1, isiz2, isiz3, nx, ny, nz,A, B, C, D, E, F)
int isiz1, isiz2, isiz3;
int nx,ny, nz;
double *A, *B, *C, *D, *E, *F;
{
  
  /* forward elimination */

  for (k=2;k<=nz-1;k++)
    {
      for (j=2;j<=ny-1;j++)
	{
	  for (i=2;i<=nx-1;i++)
	    {
	      tmp1 = b(i,j,k) / c(i-1,j,k);
	      tmp2 = a(i+1,j,k) / c(i-1,j,k);
	      c(i,j,k) = c(i,j,k) - tmp1 * d(i-1,j,k);
	      d(i,j,k) = d(i,j,k) - tmp1 * e(i-1,j,k) ;
	      f(1,i,j,k) = f(1,i,j,k) - tmp1 * f(1,i-1,j,k);
	      f(2,i,j,k) = f(2,i,j,k) - tmp1 * f(2,i-1,j,k);
	      f(3,i,j,k) = f(3,i,j,k) - tmp1 * f(3,i-1,j,k);
	      b(i+1,j,k) = b(i+1,j,k) - tmp2 * d(i-1,j,k);
	      c(i+1,j,k) = c(i+1,j,k) - tmp2 * e(i-1,j,k);
	      f(1,i+1,j,k) = f(1,i+1,j,k) - tmp2 * f(1,i-1,j,k);
	      f(2,i+1,j,k) = f(2,i+1,j,k) - tmp2 * f(2,i-1,j,k);
	      f(3,i+1,j,k) = f(3,i+1,j,k) - tmp2 * f(3,i-1,j,k);
	    }
	  tmp1 = b(nx,j,k) / c(nx-1,j,k);
	  c(nx,j,k) = c(nx,j,k) - tmp1 * d(nx-1,j,k);
	  d(nx,j,k) = d(nx,j,k) - tmp1 * e(nx-1,j,k);
	  f(1,nx,j,k) = f(1,nx,j,k) - tmp1 * f(1,nx-1,j,k);
	  f(2,nx,j,k) = f(2,nx,j,k) - tmp1 * f(2,nx-1,j,k);
	  f(3,nx,j,k) = f(3,nx,j,k) - tmp1 * f(3,nx-1,j,k);
	}
    }

  /* back-substitution phase */

  for (k=2;k<=nz-1;k++)
    {
      for (j=2;j<=ny-1;j++)
	{
	  f(1,nx,j,k) = f(1,nx,j,k) / c(nx,j,k);
	  f(1,nx-1,j,k) = ( f(1,nx-1,j,k) - d(nx-1,j,k)*f(1,nx,j,k) )
	    / c(nx-1,j,k);
	  f(2,nx,j,k) = f(2,nx,j,k) / c(nx,j,k);
	  f(2,nx-1,j,k) = ( f(2,nx-1,j,k) - d(nx-1,j,k)*f(2,nx,j,k) )
	    / c(nx-1,j,k);
	  f(3,nx,j,k) = f(3,nx,j,k) / c(nx,j,k);
	  f(3,nx-1,j,k) = ( f(3,nx-1,j,k) - d(nx-1,j,k)*f(3,nx,j,k) )
	    / c(nx-1,j,k)	;
	  for (i=nx-2;i>=1;i--)
	    {
	      f(1,i,j,k) = ( f(1,i,j,k) - d(i,j,k)*f(1,i+1,j,k)
			    - e(i,j,k)*f(1,i+2,j,k) ) / c(i,j,k);
	      f(2,i,j,k) = ( f(2,i,j,k) - d(i,j,k)*f(2,i+1,j,k)
			    - e(i,j,k)*f(2,i+2,j,k) ) / c(i,j,k);
	      f(3,i,j,k) = ( f(3,i,j,k) - d(i,j,k)*f(3,i+1,j,k)
			    - e(i,j,k)*f(3,i+2,j,k) ) / c(i,j,k);
	    }
	}
    }
}

/* solution of multiple, independent systems of penta-diagonal systems */

spentay ( isiz1, isiz2, isiz3, m, nx, ny, nz,A, B, C, D, E, F)
int isiz1, isiz2, isiz3;
int m, nx,ny, nz;
double *A, *B, *C, *D, *E, *F;
{
  
  /* forward elimination */

  for (k=2;k<=nz-1;k++)
    {
      for (i=2;i<=nx-1;i++)
	{
	  for (j=2;j<=ny-1;j++)
	    {
	      tmp1 = b(i,j,k) / c(i,j-1,k);
	      tmp2 = a(i,j+1,k) / c(i,j-1,k);
	      c(i,j,k) = c(i,j,k) - tmp1 * d(i,j-1,k);
	      d(i,j,k) = d(i,j,k) - tmp1 * e(i,j-1,k) ;
	      f(m,i,j,k) = f(m,i,j,k) - tmp1 * f(m,i,j-1,k);
	      b(i,j+1,k) = b(i,j+1,k) - tmp2 * d(i,j-1,k);
	      c(i,j+1,k) = c(i,j+1,k) - tmp2 * e(i,j-1,k);
	      f(m,i,j+1,k) = f(m,i,j+1,k) - tmp2 * f(m,i,j-1,k);
	    }
	  tmp1 = b(i,ny,k) / c(i,ny-1,k);
	  c(i,ny,k) = c(i,ny,k) - tmp1 * d(i,ny-1,k);
	  d(i,ny,k) = d(i,ny,k) - tmp1 * e(i,ny-1,k);
	  f(m,i,ny,k) = f(m,i,ny,k) - tmp1 * f(m,i,ny-1,k);
	}
    }

  /* back-substitution phase */

  for (k=2;k<=nz-1;k++)
    {
      for (i=2;i<=nx-1;i++)
	{
	  f(m,i,ny,k) = f(m,i,ny,k) / c(i,ny,k);
	  f(m,i,ny-1,k) = ( f(m,i,ny-1,k) - d(i,ny-1,k)*f(m,i,ny,k) )
	    / c(i,ny-1,k);
	  for (j=ny-2;j>=1;j--)
	    {
	      f(m,i,j,k) = ( f(m,i,j,k) - d(i,j,k)*f(m,i,j+1,k)
			    - e(i,j,k)*f(m,i,j+2,k) ) / c(i,j,k);
	    }
	}
    }
}


/* solution of multiple, independent systems of penta-diagonal systems */

spentay3 ( isiz1, isiz2, isiz3, nx, ny, nz, A, B, C, D, E, F)
int isiz1, isiz2, isiz3;
int nx,ny, nz;
double *A, *B, *C, *D, *E, *F;
{
  /* forward elimination */

  for (k=2;k<=nz-1;k++)
    {
      for (i=2;i<=nx-1;i++)
	{
	  for (j=2;j<=ny-1;j++)
	    {
	      tmp1 = b(i,j,k) / c(i,j-1,k);
	      tmp2 = a(i,j+1,k) / c(i,j-1,k);
	      c(i,j,k) = c(i,j,k) - tmp1 * d(i,j-1,k);
	      d(i,j,k) = d(i,j,k) - tmp1 * e(i,j-1,k) ;
	      f(1,i,j,k) = f(1,i,j,k) - tmp1 * f(1,i,j-1,k);
	      f(2,i,j,k) = f(2,i,j,k) - tmp1 * f(2,i,j-1,k);
	      f(3,i,j,k) = f(3,i,j,k) - tmp1 * f(3,i,j-1,k);
	      b(i,j+1,k) = b(i,j+1,k) - tmp2 * d(i,j-1,k);
	      c(i,j+1,k) = c(i,j+1,k) - tmp2 * e(i,j-1,k);
	      f(1,i,j+1,k) = f(1,i,j+1,k) - tmp2 * f(1,i,j-1,k);
	      f(2,i,j+1,k) = f(2,i,j+1,k) - tmp2 * f(2,i,j-1,k);
	      f(3,i,j+1,k) = f(3,i,j+1,k) - tmp2 * f(3,i,j-1,k);
	    }
	  tmp1 = b(i,ny,k) / c(i,ny-1,k);
	  c(i,ny,k) = c(i,ny,k) - tmp1 * d(i,ny-1,k);
	  d(i,ny,k) = d(i,ny,k) - tmp1 * e(i,ny-1,k);
	  f(1,i,ny,k) = f(1,i,ny,k) - tmp1 * f(1,i,ny-1,k);
	  f(2,i,ny,k) = f(2,i,ny,k) - tmp1 * f(2,i,ny-1,k);
	  f(3,i,ny,k) = f(3,i,ny,k) - tmp1 * f(3,i,ny-1,k);
	}
    }

  /* back-substitution phase */

  for (k=2;k<=nz-1;k++)
    {
      for (i=2;i<=nx-1;i++)
	{
	  f(1,i,ny,k) = f(1,i,ny,k) / c(i,ny,k);
	  f(1,i,ny-1,k) = ( f(1,i,ny-1,k) - d(i,ny-1,k)*f(1,i,ny,k) )
	    / c(i,ny-1,k);
	  f(2,i,ny,k) = f(2,i,ny,k) / c(i,ny,k);
	  f(2,i,ny-1,k) = ( f(2,i,ny-1,k) - d(i,ny-1,k)*f(2,i,ny,k) )
	    / c(i,ny-1,k);
	  f(3,i,ny,k) = f(3,i,ny,k) / c(i,ny,k);
	  f(3,i,ny-1,k) = ( f(3,i,ny-1,k) - d(i,ny-1,k)*f(3,i,ny,k) )
	    / c(i,ny-1,k);
	  for (j=ny-2;j>=1;j--)
	    {
	      f(1,i,j,k) = ( f(1,i,j,k) - d(i,j,k)*f(1,i,j+1,k)
			    - e(i,j,k)*f(1,i,j+2,k) ) / c(i,j,k);
	      f(2,i,j,k) = ( f(2,i,j,k) - d(i,j,k)*f(2,i,j+1,k)
			    - e(i,j,k)*f(2,i,j+2,k) ) / c(i,j,k);
	      f(3,i,j,k) = ( f(3,i,j,k) - d(i,j,k)*f(3,i,j+1,k)
			    - e(i,j,k)*f(3,i,j+2,k) ) / c(i,j,k);
	    }
	}
    }
}

/* solution of multiple, independent systems of penta-diagonal systems */

spentaz ( isiz1, isiz2, isiz3, m, nx, ny, nz, A, B, C, D, E, F)
int isiz1, isiz2, isiz3;
int m, nx,ny, nz;
double *A, *B, *C, *D, *E, *F;
{
  
  /* forward elimination */

  for (j=2;j<=ny-1;j++)
    {
      for (i=2;i<=nx-1;i++)
	{
	  for (k=2;k<=nz-1;k++)
	    {
	      tmp1 = b(i,j,k) / c(i,j,k-1);
	      tmp2 = a(i,j,k+1) / c(i,j,k-1);
	      c(i,j,k) = c(i,j,k) - tmp1 * d(i,j,k-1);
	      d(i,j,k) = d(i,j,k) - tmp1 * e(i,j,k-1) ;
	      f(m,i,j,k) = f(m,i,j,k) - tmp1 * f(m,i,j,k-1);
	      b(i,j,k+1) = b(i,j,k+1) - tmp2 * d(i,j,k-1);
	      c(i,j,k+1) = c(i,j,k+1) - tmp2 * e(i,j,k-1);
	      f(m,i,j,k+1) = f(m,i,j,k+1) - tmp2 * f(m,i,j,k-1);
	    }
	  tmp1 = b(i,j,nz) / c(i,j,nz-1);
	  c(i,j,nz) = c(i,j,nz) - tmp1 * d(i,j,nz-1);
	  d(i,j,nz) = d(i,j,nz) - tmp1 * e(i,j,nz-1);
	  f(m,i,j,nz) = f(m,i,j,nz) - tmp1 * f(m,i,j,nz-1);
	}
    }

  /* back-substitution phase */

  for (j=2;j<=ny-1;j++)
    {
      for (i=2;i<=nx-1;i++)
	{
	  f(m,i,j,nz) = f(m,i,j,nz) / c(i,j,nz);
	  f(m,i,j,nz-1) = ( f(m,i,j,nz-1) - d(i,j,nz-1)*f(m,i,j,nz) )
	    / c(i,j,nz-1);
	  for (k=nz-2;k>=1;k--)
	    {
	      f(m,i,j,k) = ( f(m,i,j,k) - d(i,j,k)*f(m,i,j,k+1)
			    - e(i,j,k)*f(m,i,j,k+2) ) / c(i,j,k);
	    }
	}
    }
}

/* solution of multiple, independent systems of penta-diagonal systems */

spentaz3 ( isiz1, isiz2, isiz3,nx, ny, nz, A, B, C, D, E, F)
int isiz1, isiz2, isiz3;
int nx,ny, nz;
double *A, *B, *C, *D, *E, *F;
{
  /* forward elimination */

  for (j=2;j<=ny-1;j++)
    {
      for (i=2;i<=nx-1;i++)
	{
	  for (k=2;k<=nz-1;k++)
	    {
	      tmp1 = b(i,j,k) / c(i,j,k-1);
	      tmp2 = a(i,j,k+1) / c(i,j,k-1);
	      c(i,j,k) = c(i,j,k) - tmp1 * d(i,j,k-1);
	      d(i,j,k) = d(i,j,k) - tmp1 * e(i,j,k-1) ;
	      f(1,i,j,k) = f(1,i,j,k) - tmp1 * f(1,i,j,k-1);
	      f(2,i,j,k) = f(2,i,j,k) - tmp1 * f(2,i,j,k-1);
	      f(3,i,j,k) = f(3,i,j,k) - tmp1 * f(3,i,j,k-1);
	      b(i,j,k+1) = b(i,j,k+1) - tmp2 * d(i,j,k-1);
	      c(i,j,k+1) = c(i,j,k+1) - tmp2 * e(i,j,k-1);
	      f(1,i,j,k+1) = f(1,i,j,k+1) - tmp2 * f(1,i,j,k-1);
	      f(2,i,j,k+1) = f(2,i,j,k+1) - tmp2 * f(2,i,j,k-1);
	      f(3,i,j,k+1) = f(3,i,j,k+1) - tmp2 * f(3,i,j,k-1);
	    }
	  tmp1 = b(i,j,nz) / c(i,j,nz-1);
	  c(i,j,nz) = c(i,j,nz) - tmp1 * d(i,j,nz-1);
	  d(i,j,nz) = d(i,j,nz) - tmp1 * e(i,j,nz-1);
	  f(1,i,j,nz) = f(1,i,j,nz) - tmp1 * f(1,i,j,nz-1);
	  f(2,i,j,nz) = f(2,i,j,nz) - tmp1 * f(2,i,j,nz-1);
	  f(3,i,j,nz) = f(3,i,j,nz) - tmp1 * f(3,i,j,nz-1);
	}
    }

  /* back-substitution phase */

  for (j=2;j<=ny-1;j++)
    {
      for (i=2;i<=nx-1;i++)
	{
	  f(1,i,j,nz) = f(1,i,j,nz) / c(i,j,nz);
	  f(1,i,j,nz-1) = ( f(1,i,j,nz-1) - d(i,j,nz-1)*f(1,i,j,nz) )
	    / c(i,j,nz-1);
	  f(2,i,j,nz) = f(2,i,j,nz) / c(i,j,nz);
	  f(2,i,j,nz-1) = ( f(2,i,j,nz-1) - d(i,j,nz-1)*f(2,i,j,nz) )
	    / c(i,j,nz-1);
	  f(3,i,j,nz) = f(3,i,j,nz) / c(i,j,nz);
	  f(3,i,j,nz-1) = ( f(3,i,j,nz-1) - d(i,j,nz-1)*f(3,i,j,nz) )
	    / c(i,j,nz-1);
	  for (k=nz-2;k>=1;k--)
	    {
	      f(1,i,j,k) = ( f(1,i,j,k) - d(i,j,k)*f(1,i,j,k+1)
			    - e(i,j,k)*f(1,i,j,k+2) ) / c(i,j,k);
	      f(2,i,j,k) = ( f(2,i,j,k) - d(i,j,k)*f(2,i,j,k+1)
			    - e(i,j,k)*f(2,i,j,k+2) ) / c(i,j,k);
	      f(3,i,j,k) = ( f(3,i,j,k) - d(i,j,k)*f(3,i,j,k+1)
			    - e(i,j,k)*f(3,i,j,k+2) ) / c(i,j,k);
	    }
	}
    }
}


