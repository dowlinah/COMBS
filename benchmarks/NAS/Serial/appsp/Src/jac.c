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

/* jac.c */

#define MAIN extern
#include "common.h"

#define cv(dim1) 	CV[dim1-1]
#define aa(dim1)	AA[dim1-1]
#define rhon(dim1)	RHON[dim1-1]
#define rhoq(dim1)	RHOQ[dim1-1]
#define rhos(dim1)	RHOS[dim1-1]
  
/* form the xi-direction pentadiagonal system */
  
jacx (m)
int m;
{
  double CV[ISIZ1];
  double AA[ISIZ1];
  double RHON[ISIZ1];
  int i, j, k;
  double r43, c34, c1345, sn, ru1, uu, vv,ww, q;
  
  r43 = 4.0 / 3.0;
  c34 = c3 * c4;
  c1345 = c1 * c3 * c4 * c5;
  if ( m == 3 )
    {
      sn = 0.0;
    }
  else if ( m == 4 )
    {
      sn = 1.0;
    }
  else if ( m == 5 )
    {
      sn = - 1.0;
    }
  
  for (k=2;k<=nz-1;k++)
    {
      for (j=2;j<=ny-1;j++)
	{
	  for (i=1;i<=nx;i++)
	    {
	      ru1 = 1.0 / u(1,i,j,k);
	      uu = ru1 * u(2,i,j,k);
	      vv = ru1 * u(3,i,j,k);
	      ww = ru1 * u(4,i,j,k);
	      q = 0.50 * ( uu*uu
			  +vv*vv
			  +ww*ww );
	      cv(i) = uu;
	      aa(i) = sqrt ( c1 * c2 * ( ru1 * u(5,i,j,k) - q ) );
	      rhon(i) = max5 ( dx1,
			      dx2 + r43 * c34 * ru1,
			      dx3 + c34 * ru1,
			      dx4 + c34 * ru1,
			      dx5 + c1345 * ru1 );
	    }
	  a(1,j,k) = 0.0;
	  b(1,j,k) = 0.0;
	  c(1,j,k) = 1.0;
	  d(1,j,k) = 0.0;
	  e(1,j,k) = 0.0;
	  for (i=2;i<=nx-1;i++)
	    {
	      a(i,j,k) =   0.0;
	      b(i,j,k) = - dt * tx2 * ( cv(i-1) + sn * aa(i-1) )
		- dt * rhon(i-1) * tx1 ;
	      c(i,j,k) =   1.0
		+ dt * rhon(i) * tx1 * 2.0;
	      d(i,j,k) =   dt * tx2 * ( cv(i+1) + sn * aa(i+1) )
		- dt * rhon(i+1) * tx1;
	      e(i,j,k) =   0.0;
	    }
	  a(nx,j,k) = 0.0;
	  b(nx,j,k) = 0.0;
	  c(nx,j,k) = 1.0;
	  d(nx,j,k) = 0.0;
	  e(nx,j,k) = 0.0;
	  
	  /* fourth order dissipation  */
	  c(2,j,k) = c(2,j,k) + dt * dssp * (  5.0 );
	  d(2,j,k) = d(2,j,k) + dt * dssp * ( -4.0 );
	  e(2,j,k) = e(2,j,k) + dt * dssp * (  1.0 );
	  b(3,j,k) = b(3,j,k) + dt * dssp * ( -4.0 );
	  c(3,j,k) = c(3,j,k) + dt * dssp * ( 6.0 );
	  d(3,j,k) = d(3,j,k) + dt * dssp * ( -4.0 );
	  e(3,j,k) = e(3,j,k) + dt * dssp * ( 1.0 );
	  for (i=4;i<=nx-3;i++)
	    {
	      a(i,j,k) = a(i,j,k) + dt * dssp * ( 1.0 );
	      b(i,j,k) = b(i,j,k) + dt * dssp * ( -4.0 );
	      c(i,j,k) = c(i,j,k) + dt * dssp * ( 6.0 );
	      d(i,j,k) = d(i,j,k) + dt * dssp * ( -4.0 );
	      e(i,j,k) = e(i,j,k) + dt * dssp * ( 1.0 );
	    }
	  a(nx-2,j,k) = a(nx-2,j,k) + dt * dssp * ( 1.0 );
	  b(nx-2,j,k) = b(nx-2,j,k) + dt * dssp * ( -4.0 );
	  c(nx-2,j,k) = c(nx-2,j,k) + dt * dssp * ( 6.0 );
	  d(nx-2,j,k) = d(nx-2,j,k) + dt * dssp * ( -4.0 );
	  a(nx-1,j,k) = a(nx-1,j,k) + dt * dssp * ( 1.0 );
	  b(nx-1,j,k) = b(nx-1,j,k) + dt * dssp * ( -4.0 );
	  c(nx-1,j,k) = c(nx-1,j,k) + dt * dssp * ( 5.0 );
	}
    }
}

/* form the eta-direction pentadiagonal system. */

jacy(m)
int m;
{
  double CV[ISIZ2];
  double AA[ISIZ2];
  double RHOQ[ISIZ2];
  int i, j, k;
  double r43, c34, c1345, sn, ru1, uu, vv, ww, q;
  
  r43 = 4.0 / 3.0;
  c34 = c3 * c4;
  c1345 = c1 * c3 * c4 * c5;
  if ( m == 3 )
    {
      sn = 0.0;
    }
  else if ( m == 4 )
    {
      sn = 1.0;
    }
  else if ( m == 5 )
    {
      sn = -1.0;
    }   
  for (k=2;k<=nz-1;k++)
    {
      for (i=2;i<=nx-1;i++)
	{
	  for (j=1;j<=ny;j++)
	    {
	      ru1 = 1.0 / u(1,i,j,k);
	      uu = ru1 * u(2,i,j,k);
	      vv = ru1 * u(3,i,j,k);
	      ww = ru1 * u(4,i,j,k);
	      q = 0.50 * ( uu*uu
			  +vv*vv
			  +ww*ww );
	      cv(j) = vv;
	      aa(j) = sqrt ( c1 * c2 * ( ru1 * u(5,i,j,k) - q ) );
	      rhoq(j) = max5 ( dy1,
			      dy2 + c34 * ru1,
			      dy3 + r43 * c34 * ru1,
			      dy4 + c34 * ru1,
			      dy5 + c1345 * ru1 );
	    }
	  a(i,1,k) = 0.0;
	  b(i,1,k) = 0.0;
	  c(i,1,k) = 1.0;
	  d(i,1,k) = 0.0;
	  e(i,1,k) = 0.0;
	  for (j=2;j<=ny-1;j++)
	    {
	      a(i,j,k) = 0.0;
	      b(i,j,k) = - dt * ty2 * ( cv(j-1) + sn * aa(j-1) )
		- dt * rhoq(j-1) * ty1 ;
	      c(i,j,k) =   1.0
		+ dt * rhoq(j) * ty1 * 2.0;
	      d(i,j,k) =   dt * ty2 * ( cv(j+1) + sn * aa(j+1) )
		- dt * rhoq(j+1) * ty1;
	      e(i,j,k) = 0.0;
	    }
	  a(i,ny,k) = 0.0;
	  b(i,ny,k) = 0.0;
	  c(i,ny,k) = 1.0;
	  d(i,ny,k) = 0.0;
	  e(i,ny,k) = 0.0;
	  
	  /* fourth order dissipation  */
	  
	  c(i,2,k) = c(i,2,k) + dt * dssp * ( 5.0 );
	  d(i,2,k) = d(i,2,k) + dt * dssp * ( -4.0 );
	  e(i,2,k) = e(i,2,k) + dt * dssp * ( 1.0 );
	  b(i,3,k) = b(i,3,k) + dt * dssp * ( -4.0 );
	  c(i,3,k) = c(i,3,k) + dt * dssp * ( 6.0 );
	  d(i,3,k) = d(i,3,k) + dt * dssp * ( -4.0 );
	  e(i,3,k) = e(i,3,k) + dt * dssp * ( 1.0 );
	  for (j=4;j<=ny-3;j++)
	    {
	      a(i,j,k) = a(i,j,k) + dt * dssp * ( 1.0 );
	      b(i,j,k) = b(i,j,k) + dt * dssp * ( -4.0 );
	      c(i,j,k) = c(i,j,k) + dt * dssp * ( 6.0 );
	      d(i,j,k) = d(i,j,k) + dt * dssp * ( -4.0 );
	      e(i,j,k) = e(i,j,k) + dt * dssp * ( 1.0 );
	    }
	  a(i,ny-2,k) = a(i,ny-2,k) + dt * dssp * ( 1.0 );
	  b(i,ny-2,k) = b(i,ny-2,k) + dt * dssp * ( -4.0 );
	  c(i,ny-2,k) = c(i,ny-2,k) + dt * dssp * ( 6.0 );
	  d(i,ny-2,k) = d(i,ny-2,k) + dt * dssp * ( -4.0 );
	  a(i,ny-1,k) = a(i,ny-1,k) + dt * dssp * ( 1.0 );
	  b(i,ny-1,k) = b(i,ny-1,k) + dt * dssp * ( -4.0 );
	  c(i,ny-1,k) = c(i,ny-1,k) + dt * dssp * ( 5.0 );
	}
    }
}

/* form the zeta-direction pentadiagonal system. */

jacz(m)
int m;
{
  double CV[ISIZ3];
  double AA[ISIZ3];
  double RHOS[ISIZ3];
  int i, j, k;
  double r43, c34, c1345, sn, ru1, uu, vv, ww, q;  
  
  r43 = 4.0 / 3.0;
  c34 = c3 * c4;
  c1345 = c1 * c3 * c4 * c5;
  if ( m == 3 )
    {
      sn = 0.0;
    }
  else if ( m == 4 )
    {
      sn = 1.0;
    }
  else if ( m == 5 )
    {
      sn = - 1.0;
    }   
  for (j=2;j<=ny-1;j++)
    {
      for (i=2;i<=nx-1;i++)
	{
	  for (k=1;k<=nz;k++)
	    {
	      ru1 = 1.0 / u(1,i,j,k);
	      uu = ru1 * u(2,i,j,k);
	      vv = ru1 * u(3,i,j,k);
	      ww = ru1 * u(4,i,j,k);
	      q = 0.50 * ( uu*uu
			  +vv*vv
			  +ww*ww );
	      cv(k) = ww;
	      aa(k) = sqrt ( c1 * c2 * ( ru1 * u(5,i,j,k) - q ) );
	      rhos(k) = max5 ( dz1,
			      dz2 + c34 * ru1 ,
			      dz3 + c34 * ru1 ,
			      dz4 + r43 * c34 * ru1,
			      dz5 + c1345 * ru1 );
	    }
	  a(i,j,1) = 0.0;
	  b(i,j,1) = 0.0;
	  c(i,j,1) = 1.0;
	  d(i,j,1) = 0.0;
	  e(i,j,1) = 0.0;
	  for (k=2;k<=nz-1;k++)
	    {
	      a(i,j,k) = 0.0;
	      b(i,j,k) = - dt * tz2 * ( cv(k-1) + sn * aa(k-1) )
		- dt * rhos(k-1) * tz1 ;
	      c(i,j,k) =   1.0
		+ dt * rhos(k) * tz1 * 2.0;
	      d(i,j,k) =   dt * tz2 * ( cv(k+1) + sn * aa(k+1) )
		- dt * rhos(k+1) * tz1;
	      e(i,j,k) = 0.0;
	    }
	  a(i,j,nz) = 0.0;
	  b(i,j,nz) = 0.0;
	  c(i,j,nz) = 1.0;
	  d(i,j,nz) = 0.0;
	  e(i,j,nz) = 0.0;
	  
	  /* fourth order dissipation  */
	  
	  c(i,j,2) = c(i,j,2) + dt * dssp * ( 5.0 );
	  d(i,j,2) = d(i,j,2) + dt * dssp * ( -4.0 );
	  e(i,j,2) = e(i,j,2) + dt * dssp * ( 1.0 );
	  b(i,j,3) = b(i,j,3) + dt * dssp * ( -4.0 );
	  c(i,j,3) = c(i,j,3) + dt * dssp * ( 6.0 );
	  d(i,j,3) = d(i,j,3) + dt * dssp * ( -4.0 );
	  e(i,j,3) = e(i,j,3) + dt * dssp * ( 1.0 );
	  for (k=4;k<=nz-3;k++)
	    {
	      a(i,j,k) = a(i,j,k) + dt * dssp * ( 1.0 );
	      b(i,j,k) = b(i,j,k) + dt * dssp * ( -4.0 );
	      c(i,j,k) = c(i,j,k) + dt * dssp * ( 6.0 );
	      d(i,j,k) = d(i,j,k) + dt * dssp * ( -4.0 );
	      e(i,j,k) = e(i,j,k) + dt * dssp * ( 1.0 );
	    }
	  a(i,j,nz-2) = a(i,j,nz-2) + dt * dssp * ( 1.0 );
	  b(i,j,nz-2) = b(i,j,nz-2) + dt * dssp * ( -4.0 );
	  c(i,j,nz-2) = c(i,j,nz-2) + dt * dssp * ( 6.0 );
	  d(i,j,nz-2) = d(i,j,nz-2) + dt * dssp * ( -4.0 );
	  a(i,j,nz-1) = a(i,j,nz-1) + dt * dssp * (  1.0 );
	  b(i,j,nz-1) = b(i,j,nz-1) + dt * dssp * ( -4.0 );
	  c(i,j,nz-1) = c(i,j,nz-1) + dt * dssp * ( 5.0 );
	}
    }
}
