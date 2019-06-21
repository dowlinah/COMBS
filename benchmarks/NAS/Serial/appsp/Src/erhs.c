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

/* erhs.c */

#define MAIN extern
#include "common.h"


#define flux(dim1, dim2)	TWO_D_ARRAY(FLUX, dim1, dim2)
#define ue(dim1, dim2)		TWO_D_ARRAY(UE, dim1, dim2)

/* compute the right hand side based on exact solution */

erhs()
{
  double FLUX[TWO_D_SIZE], UE[TWO_D_SIZE];
  double tmp;
  double dsspm = dssp;
  int i,j,k,m;
  double xi, eta, zeta;
  double q, u21, u31, u41, u51, u21i, u31i, u41i, u51i;
  double u21im1, u31im1, u41im1, u51im1;
  double u21j, u31j, u41j, u51j;
  double u21jm1, u31jm1, u41jm1, u51jm1;
  double u21k, u31k, u41k, u51k;
  double u21km1, u31km1, u41km1, u51km1;

  for (k=1;k<=nz;k++) 
    {
      for (j=1;j<=ny;j++) 
	{
	  for (i=1;i<=nx;i++) 
	    {
	      for (m=1;m<=5;m++)
		{
		  frct( m, i, j, k ) = 0.0;
		}
            }
	}
    }
  
  /* xi-direction flux differences */

  for (k=2;k<=nz-1;k++) 
    {
      zeta = ( (double)(k-1) ) / ( nz - 1 );
      for (j=2;j<=ny-1;j++) 
	{
	  eta = ( (double)(j-1) ) / ( ny - 1 );
	  for (i=1;i<=nx;i++) 
	    {
	      xi = ( (double)(i-1) ) / ( nx - 1 );
	      for (m=1;m<=5;m++) 
		{
                  ue(m,i) =  ce(m,1)
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
	      flux(1,i) = ue(2,i);
	      u21 = ue(2,i) / ue(1,i);
	      q = 0.50 * (  ue(2,i) * ue(2,i)
			  + ue(3,i) * ue(3,i)
			  + ue(4,i) * ue(4,i) )
		/ ue(1,i);
	      flux(2,i) = ue(2,i) * u21 + c2 * ( ue(5,i) - q );
	      flux(3,i) = ue(3,i) * u21;
	      flux(4,i) = ue(4,i) * u21;
	      flux(5,i) = ( c1 * ue(5,i) - c2 * q ) * u21;
            }
	  for (i=2;i<=nx-1;i++) 
	    {
	      for (m=1;m<=5;m++) 
		{
                  frct(m,i,j,k) =  frct(m,i,j,k)
		    - tx2 * ( flux(m,i+1) - flux(m,i-1) );
		}
            }
	  for (i=2;i<=nx;i++) 
	    {
	      tmp = 1.0 / ue(1,i);
	      u21i = tmp * ue(2,i);
	      u31i = tmp * ue(3,i);
	      u41i = tmp * ue(4,i);
	      u51i = tmp * ue(5,i);
	      tmp = 1.0 / ue(1,i-1);
	      u21im1 = tmp * ue(2,i-1);
	      u31im1 = tmp * ue(3,i-1);
	      u41im1 = tmp * ue(4,i-1);
	      u51im1 = tmp * ue(5,i-1);
	      flux(2,i) = (4.0/3.0) * tx3 * ( u21i - u21im1 );
	      flux(3,i) = tx3 * ( u31i - u31im1 );
	      flux(4,i) = tx3 * ( u41i - u41im1 );
	      flux(5,i) = 0.50 * ( 1.0 - c1*c5 )
		* tx3 * ( ( u21i*u21i + u31i*u31i + u41i*u41i )
			 - (u21im1*u21im1 +u31im1*u31im1 +u41im1*u41im1 ) )
		  + (1.0/6.0)
		    * tx3 * (u21i*u21i -u21im1*u21im1 )
		      + c1 * c5 * tx3 * ( u51i - u51im1 );
            }

	  for (i=2;i<=nx-1;i++) 
	    {
	      frct(1,i,j,k) = frct(1,i,j,k)
		+ dx1 * tx1 * (            ue(1,i-1)
			       - 2.0 * ue(1,i)
			       +           ue(1,i+1) );
	      frct(2,i,j,k) = frct(2,i,j,k)
		+ tx3 * c3 * c4 * ( flux(2,i+1) - flux(2,i) )
		  + dx2 * tx1 * (            ue(2,i-1)
				 - 2.0 * ue(2,i)
				 +           ue(2,i+1) );
	      frct(3,i,j,k) = frct(3,i,j,k)
		+ tx3 * c3 * c4 * ( flux(3,i+1) - flux(3,i) )
		  + dx3 * tx1 * (            ue(3,i-1)
				 - 2.0 * ue(3,i)
				 +           ue(3,i+1) );
	      frct(4,i,j,k) = frct(4,i,j,k)
		+ tx3 * c3 * c4 * ( flux(4,i+1) - flux(4,i) )
		  + dx4 * tx1 * (            ue(4,i-1)
                                  - 2.0 * ue(4,i)
                                  +           ue(4,i+1) );
	      frct(5,i,j,k) = frct(5,i,j,k)
		+ tx3 * c3 * c4 * ( flux(5,i+1) - flux(5,i) )
		  + dx5 * tx1 * (            ue(5,i-1)
                                  - 2.0 * ue(5,i)
                                   +           ue(5,i+1) );
            }
	  
	  /* Fourth-order dissipation */

	  for (m=1;m<=5;m++) 
	    {
	      frct(m,2,j,k) = frct(m,2,j,k)
                - dsspm * (  5.0 * ue(m,2)
			   - 4.0 * ue(m,3)
			   +           ue(m,4) );
	      frct(m,3,j,k) = frct(m,3,j,k)
                - dsspm * ( - 4.0 * ue(m,2)
			   + 6.0 * ue(m,3)
			   - 4.0 * ue(m,4)
			   +           ue(m,5) );
            }
	  for (i=4;i<=nx-3;i++) 
	    {
	      for (m=1;m<=5;m++) 
		{
                  frct(m,i,j,k) = frct(m,i,j,k)
		    - dsspm * (            ue(m,i-2)
			       - 4.0 * ue(m,i-1)
			       + 6.0 * ue(m,i)
			       - 4.0 * ue(m,i+1)
			       +           ue(m,i+2) );
		}
            }
	  for (m=1;m<=5;m++) 
	    {
	      frct(m,nx-2,j,k) = frct(m,nx-2,j,k)
                - dsspm * (             ue(m,nx-4)
			   - 4.0 * ue(m,nx-3)
			   + 6.0 * ue(m,nx-2)
			   - 4.0 * ue(m,nx-1)  );
	      frct(m,nx-1,j,k) = frct(m,nx-1,j,k)
                - dsspm * (             ue(m,nx-3)
			   - 4.0 * ue(m,nx-2)
			   + 5.0 * ue(m,nx-1) );
            }
	}
    }
  
  /* eta-direction flux differences */

  for (k=2;k<=nz-1;k++) 
    {
      zeta = ( (double)(k-1) ) / ( nz - 1 );
      for (i=2;i<=nx-1;i++) 	
	{
	  xi = ( (double)(i-1) ) / ( nx - 1 );
	  for (j=1;j<=ny;j++) {
	    eta = ( (double)(j-1) ) / ( ny - 1 );
	    for (m=1;m<=5;m++) 
	      {
		ue(m,j) =  ce(m,1)
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
	    flux(1,j) = ue(3,j);
	    u31 = ue(3,j) / ue(1,j);
	    q = 0.50 * (  ue(2,j) * ue(2,j)
			+ ue(3,j) * ue(3,j)
			+ ue(4,j) * ue(4,j) )
	      / ue(1,j);
	    flux(2,j) = ue(2,j) * u31;
	    flux(3,j) = ue(3,j) * u31 + c2 * ( ue(5,j) - q );
	    flux(4,j) = ue(4,j) * u31;
	    flux(5,j) = ( c1 * ue(5,j) - c2 * q ) * u31;
	  }
	  for (j=2;j<=ny-1;j++) 
	    {
	      for (m=1;m<=5;m++) {
		frct(m,i,j,k) =  frct(m,i,j,k)
		  - ty2 * ( flux(m,j+1) - flux(m,j-1) );
	      }
	    }
	  for (j=2;j<=ny;j++) 
	    {
	      tmp = 1.0 / ue(1,j);
	      u21j = tmp * ue(2,j);
	      u31j = tmp * ue(3,j);
	      u41j = tmp * ue(4,j);
	      u51j = tmp * ue(5,j);
	      tmp = 1.0 / ue(1,j-1);
	      u21jm1 = tmp * ue(2,j-1);
	      u31jm1 = tmp * ue(3,j-1);
	      u41jm1 = tmp * ue(4,j-1);
	      u51jm1 = tmp * ue(5,j-1);
	      flux(2,j) = ty3 * ( u21j - u21jm1 );
	      flux(3,j) = (4.0/3.0) * ty3 * ( u31j - u31jm1 );
	      flux(4,j) = ty3 * ( u41j - u41jm1 );
	      flux(5,j) = 0.50 * ( 1.0 - c1*c5 )
		* ty3 * ( ( u21j*u21j + u31j*u31j + u41j*u41j )
			 - (u21jm1*u21jm1 +u31jm1*u31jm1 +u41jm1*u41jm1 ) )
		  + (1.0/6.0)
		    * ty3 * (u31j*u31j -u31jm1*u31jm1 )
		      + c1 * c5 * ty3 * ( u51j - u51jm1 );
            }
	  for (j=2;j<=ny-1;j++) {
	    frct(1,i,j,k) = frct(1,i,j,k)
	      + dy1 * ty1 * (            ue(1,j-1)
			     - 2.0 * ue(1,j)
			     +           ue(1,j+1) );
	    frct(2,i,j,k) = frct(2,i,j,k)
	      + ty3 * c3 * c4 * ( flux(2,j+1) - flux(2,j) )
		+ dy2 * ty1 * (            ue(2,j-1)
			       - 2.0 * ue(2,j)
			       +           ue(2,j+1) );
	    frct(3,i,j,k) = frct(3,i,j,k)
	      + ty3 * c3 * c4 * ( flux(3,j+1) - flux(3,j) )
		+ dy3 * ty1 * (            ue(3,j-1)
			       - 2.0 * ue(3,j)
			       +           ue(3,j+1) );
	    frct(4,i,j,k) = frct(4,i,j,k)
	      + ty3 * c3 * c4 * ( flux(4,j+1) - flux(4,j) )
		+ dy4 * ty1 * (            ue(4,j-1)
			       - 2.0 * ue(4,j)
			       +           ue(4,j+1) );
	    frct(5,i,j,k) = frct(5,i,j,k)
	      + ty3 * c3 * c4 * ( flux(5,j+1) - flux(5,j) )
		+ dy5 * ty1 * (            ue(5,j-1)
			       - 2.0 * ue(5,j)
			       +           ue(5,j+1) );
	  }

	  /* fourth-order dissipation */

            for (m=1;m<=5;m++) {
               frct(m,i,2,k) = frct(m,i,2,k)
                - dsspm * (  5.0 * ue(m,2)
                            - 4.0 * ue(m,3)
                            +           ue(m,4) );
               frct(m,i,3,k) = frct(m,i,3,k)
                - dsspm * ( - 4.0 * ue(m,2)
                            + 6.0 * ue(m,3)
                            - 4.0 * ue(m,4)
                            +           ue(m,5) );
            }
	    for (j=4;j<=ny-3;j++) {
               for (m=1;m<=5;m++) {
                  frct(m,i,j,k) = frct(m,i,j,k)
                   - dsspm * (            ue(m,j-2)
                             - 4.0 * ue(m,j-1)
                             + 6.0 * ue(m,j)
                             - 4.0 * ue(m,j+1)
                             +           ue(m,j+2) );
               }
            }
            for (m=1;m<=5;m++) 
	      {
		frct(m,i,ny-2,k) = frct(m,i,ny-2,k)
		  - dsspm * (             ue(m,ny-4)
			     - 4.0 * ue(m,ny-3)
			     + 6.0 * ue(m,ny-2)
			     - 4.0 * ue(m,ny-1)  );
		frct(m,i,ny-1,k) = frct(m,i,ny-1,k)
		  - dsspm * (             ue(m,ny-3)
			     - 4.0 * ue(m,ny-2)
			     + 5.0 * ue(m,ny-1)  );
	      }
	}
    }

  /* zeta-direction flux differences */
  for (j=2;j<=ny-1;j++) 
    {
      eta = ( (double)(j-1) ) / ( ny - 1 );
      for (i=2;i<=nx-1;i++) 
	{
	  xi = ( (double)(i-1) ) / ( nx - 1 );
	  for (k=1;k<=nz;k++) 
	    {
	      zeta = ( (double)(k-1) ) / ( nz - 1 );
	      for (m=1;m<=5;m++) 
		{
                  ue(m,k) =  ce(m,1)
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
	      flux(1,k) = ue(4,k);
	      u41 = ue(4,k) / ue(1,k);
	      q = 0.50 * (  ue(2,k) * ue(2,k)
			  + ue(3,k) * ue(3,k)
			  + ue(4,k) * ue(4,k) )
		/ ue(1,k);
	      flux(2,k) = ue(2,k) * u41 ;
	      flux(3,k) = ue(3,k) * u41 ;
	      flux(4,k) = ue(4,k) * u41 + c2 * ( ue(5,k) - q );
	      flux(5,k) = ( c1 * ue(5,k) - c2 * q ) * u41;
            }
	  for (k=2;k<=nz-1;k++) 
	    {
	      for (m=1;m<=5;m++) {
		frct(m,i,j,k) =  frct(m,i,j,k)
		  - tz2 * ( flux(m,k+1) - flux(m,k-1) );
	      }
            }
	  for (k=2;k<=nz;k++) 
	    {
	      tmp = 1.0 / ue(1,k);
	      u21k = tmp * ue(2,k);
	      u31k = tmp * ue(3,k);
	      u41k = tmp * ue(4,k);
	      u51k = tmp * ue(5,k);
	      tmp = 1.0 / ue(1,k-1);
	      u21km1 = tmp * ue(2,k-1);
	      u31km1 = tmp * ue(3,k-1);
	      u41km1 = tmp * ue(4,k-1);
	      u51km1 = tmp * ue(5,k-1);
	      flux(2,k) = tz3 * ( u21k - u21km1 );
	      flux(3,k) = tz3 * ( u31k - u31km1 );
	      flux(4,k) = (4.0/3.0) * tz3 * ( u41k - u41km1 );
	      flux(5,k) = 0.50 * ( 1.0 - c1*c5 )
		* tz3 * ( ( u21k*u21k + u31k*u31k + u41k*u41k )
			 - (u21km1*u21km1 +u31km1*u31km1 +u41km1*u41km1 ) )
		  + (1.0 /6.0)
		    * tz3 * (u41k*u41k -u41km1*u41km1 )
		      + c1 * c5 * tz3 * ( u51k - u51km1 );
            }
            for (k=2;k<=nz-1;k++) 
	      {
		frct(1,i,j,k) = frct(1,i,j,k)
		  + dz1 * tz1 * (            ue(1,k+1)
				 - 2.0 * ue(1,k)
				 +           ue(1,k-1) );
		frct(2,i,j,k) = frct(2,i,j,k)
		  + tz3 * c3 * c4 * ( flux(2,k+1) - flux(2,k) )
		    + dz2 * tz1 * (            ue(2,k+1)
				   - 2.0 * ue(2,k)
				   +           ue(2,k-1) );
		frct(3,i,j,k) = frct(3,i,j,k)
		  + tz3 * c3 * c4 * ( flux(3,k+1) - flux(3,k) )
		    + dz3 * tz1 * (            ue(3,k+1)
				   - 2.0 * ue(3,k)
				   +           ue(3,k-1) );
		frct(4,i,j,k) = frct(4,i,j,k)
		  + tz3 * c3 * c4 * ( flux(4,k+1) - flux(4,k) )
		    + dz4 * tz1 * (            ue(4,k+1)
				   - 2.0 * ue(4,k)
				   +           ue(4,k-1) );
		frct(5,i,j,k) = frct(5,i,j,k)
		  + tz3 * c3 * c4 * ( flux(5,k+1) - flux(5,k) )
		    + dz5 * tz1 * (            ue(5,k+1)
				   - 2.0 * ue(5,k)
				   +           ue(5,k-1) );
	      }

	  /* fourth-order dissipation */

	  for (m=1;m<=5;m++) 
	    {
	      frct(m,i,j,2) = frct(m,i,j,2)
                - dsspm * (  5.0 * ue(m,2)
			   - 4.0 * ue(m,3)
			   +           ue(m,4) );
	      frct(m,i,j,3) = frct(m,i,j,3)
                - dsspm * (- 4.0 * ue(m,2)
                           + 6.0 * ue(m,3)
                           - 4.0 * ue(m,4)
                           +           ue(m,5) );
            }
	  for (k=4;k<=nz-3;k++) 
	    {
               for (m=1;m<=5;m++) 
		 {
		   frct(m,i,j,k) = frct(m,i,j,k)
		     - dsspm * (           ue(m,k-2)
				- 4.0 * ue(m,k-1)
				+ 6.0 * ue(m,k)
				- 4.0 * ue(m,k+1)
				+           ue(m,k+2) );
		 }
	     }
	  for (m=1;m<=5;m++) 
	    {
	      frct(m,i,j,nz-2) = frct(m,i,j,nz-2)
                - dsspm * (            ue(m,nz-4)
                           - 4.0 * ue(m,nz-3)
                           + 6.0 * ue(m,nz-2)
                           - 4.0 * ue(m,nz-1)  );
	      frct(m,i,j,nz-1) = frct(m,i,j,nz-1)
                - dsspm * (             ue(m,nz-3)
			   - 4.0 * ue(m,nz-2)
			   + 5.0 * ue(m,nz-1)  );
            }
	}
    }
}




