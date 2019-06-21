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

/* adi.c */

#define MAIN extern
#include "common.h"

#define one 1.0
#define idmax(dim1) 	ONE_D_ARRAY(IDMAX, dim1)
#define jdmax(dim1)	ONE_D_ARRAY(JDMAX, dim1)
#define kdmax(dim1)	ONE_D_ARRAY(KDMAX, dim1)
#define imax(dim1)	ONE_D_ARRAY(IMAX, dim1)
#define jmax(dim1) 	ONE_D_ARRAY(JMAX, dim1)
#define kmax(dim1) 	ONE_D_ARRAY(KMAX, dim1)
#define delunm(dim1)	ONE_D_ARRAY(DELUNM, dim1)

#define v(dim1, dim2, dim3, dim4) 	FOUR_D_ARRAY(V, dim1, dim2, dim3, dim4)
#define vnm(dim1)	ONE_D_ARRAY(VNM, dim1)

/* perform pseudo-time stepping iterations for five coupled, nonlinear pde's */

adi()
{
  int IDMAX[5], JDMAX[5], KDMAX[5], IMAX[5], JMAX[5], KMAX[5];
  double DELUNM[5];
  int i,j,k,m;
  double tstart;
  int istep;
  double tend;
  double t1;
  int lnorm = 2;
  int wwt_i;

  /* compute the steady-state residuals */
   
  rhs();

  /*compute the norms of the residuals */

  if ( lnorm == 1 ) 
    {
      maxnorm ( isiz1, isiz2, isiz3,
	       nx, ny, nz,
	       &imax(1), &jmax(1), &kmax(1),
	       &rsd(1,1,1,1), &rsdnm(1) );

      if ( ipr == 1 )  
	{
	  fprintf(outfp, "          Initial Residual norms\n");
	  fprintf(outfp, "\n");

	  fprintf(outfp, " max-norm of steady-state residual ");
	  fprintf(outfp, "for first pde  = (%4f,%4f,%4f)\n", 
		  imax(1), jmax(1), kmax(1));
	  fprintf(outfp, " max-norm of steady-state residual ");
	  fprintf(outfp, "for second pde = (%4f,%4f,%4f)\n", 
		  imax(2), jmax(2), kmax(2));
	  fprintf(outfp, " max-norm of steady-state residual ");
	  fprintf(outfp, "for third pde  = (%4f,%4f,%4f)\n", 
		  imax(3), jmax(3), kmax(3));
	  fprintf(outfp, " max-norm of steady-state residual ");
	  fprintf(outfp, "for fourth pde = (%4f,%4f,%4f)\n", 
		  imax(4), jmax(4), kmax(4));
	  fprintf(outfp, " max-norm of steady-state residual ");
	  fprintf(outfp, "for fifth pde  = (%4f,%4f,%4f)\n", 
		  imax(5), jmax(5), kmax(5));
	  fflush(outfp);
	}
    }
  else if ( lnorm == 2 ) 
    {
      l2norm ( isiz1, isiz2, isiz3,
	      nx, ny, nz,
	      &rsd(1,1,1,1), &rsdnm(1) );
      if ( ipr == 1 ) 
	{
	  fprintf(outfp, "          Initial Residual norms\n");
	  fprintf(outfp, "\n");

	  fprintf(outfp, " RMS-norm of steady-state residual for ");
	  fprintf(outfp, "first pde   = %4e\n", rsdnm(1));
	  fprintf(outfp, " RMS-norm of steady-state residual for ");
	  fprintf(outfp, "second pde  = %4e\n", rsdnm(2));
	  fprintf(outfp, " RMS-norm of steady-state residual for ");
	  fprintf(outfp, "third pde   = %4e\n", rsdnm(3));
	  fprintf(outfp, " RMS-norm of steady-state residual for ");
	  fprintf(outfp, "fourth pde  = %4e\n", rsdnm(4));
	  fprintf(outfp, " RMS-norm of steady-state residual for ");
	  fprintf(outfp, "fifth pde   = %4e\n\n\n", rsdnm(5));
	  fflush(outfp);
	}
    }

  /*begin pseudo-time stepping iterations */

  tstart = timer();
  
  for (istep =1; istep<=itmax;istep++)
    {
      if ( (  mod( istep, inorm ) == 0 ) && ( ipr == 1 ) )
	{
	  fprintf(outfp, "     pseudo-time Scalar ADI iteration no. = %4d\n\n\n", istep);
	}
      for (k=1; k<=nz;k++) 
	{
	  for (j=1; j<=ny;j++) 
	    {
	      for (i=1; i<=nx;i++) 
		{
		  for (m=1;m<=5;m++) 
		    {
		      rsd(m,i,j,k) = dt * rsd(m,i,j,k);
		    }
		}
            }
	}
      
      /*perform 3-factor, scalar ADI iterations */

      /*perform the block diagonal inversion*/
      
      txinvr();

      /*perform the xsi-direction sweep */
      
      jacx ( 3 );
      spentax3 ( isiz1, isiz2, isiz3, nx, ny, nz,
		&a(1, 1, 1), &b(1, 1, 1), &c(1, 1, 1), 
		&d(1, 1, 1), &e(1, 1, 1), &rsd(1,1,1,1) );

      jacx ( 4 );
      spentax ( isiz1, isiz2, isiz3,
	       4, nx, ny, nz,
	       &a(1,1,1), &b(1,1,1), &c(1,1,1), 
	       &d(1,1,1), &e(1,1,1), &rsd(1,1,1,1) );

      jacx ( 5 );
      spentax ( isiz1, isiz2, isiz3,
	       5, nx, ny, nz,
	       &a(1,1,1), &b(1,1,1), &c(1,1,1), 
	       &d(1,1,1), &e(1,1,1), &rsd(1,1,1,1));

      /*perform the block diagonal inversion */
      ninvr();

      /*perform the eta-direction sweep */
      
      jacy ( 3 ) ;
      spentay3 ( isiz1, isiz2, isiz3,
		nx, ny, nz,
		&a(1,1,1), &b(1,1,1), &c(1,1,1), 
		&d(1,1,1), &e(1,1,1), &rsd(1,1,1,1));

      jacy ( 4 ) ;
      spentay ( isiz1, isiz2, isiz3,
	       4, nx, ny, nz,
	       &a(1,1,1), &b(1,1,1), &c(1,1,1), 
	       &d(1,1,1), &e(1,1,1), &rsd(1,1,1,1));

      jacy ( 5 ) ;
      spentay ( isiz1, isiz2, isiz3,
	       5, nx, ny, nz,
	       &a(1,1,1), &b(1,1,1), &c(1,1,1), 
	       &d(1,1,1), &e(1,1,1), &rsd(1,1,1,1));

      /*perform the block diagonal inversion */
      
      pinvr();

      /*perform the zeta-direction sweep */
      
      jacz ( 3 );
      spentaz3 ( isiz1, isiz2, isiz3,
		nx, ny, nz,
		&a(1,1,1), &b(1,1,1), &c(1,1,1), 
		&d(1,1,1), &e(1,1,1), &rsd(1,1,1,1));

      jacz ( 4 );
      spentaz ( isiz1, isiz2, isiz3,
	       4, nx, ny, nz,
	       &a(1,1,1), &b(1,1,1), &c(1,1,1), 
	       &d(1,1,1), &e(1,1,1), &rsd(1,1,1,1));

      jacz ( 5 );
      spentaz ( isiz1, isiz2, isiz3,
	       5, nx, ny, nz,
	       &a(1,1,1), &b(1,1,1), &c(1,1,1), 
	       &d(1,1,1), &e(1,1,1), &rsd(1,1,1,1));

      /*perform the block diagonal matrix-vector multiplication */
      tzetar();

      /*update the variables */

      for (k=2;k<=nz-1;k++) 
	{
	  for (j=2;j<=ny-1;j++) 
	    {
	      for (i=2;i<=nx-1;i++) 
		{
		  for (m=1;m<=5;m++) 
		    {
		      u( m, i, j, k ) =  u( m, i, j, k )
			                 + rsd( m, i, j, k );
		    }
		}
            }
	}

      /* compute the norms of pseudo-time iteration corrections */
      if ( mod ( istep, inorm ) == 0 ) 
	{
	  if ( lnorm == 1 ) 
	    {
	      maxnorm ( isiz1, isiz2, isiz3,
		       nx, ny, nz,
		       &idmax(1), &jdmax(1), &kdmax(1),
		       &rsd(1,1,1,1), &delunm(1) );
	      if ( ipr == 1 ) 
		{
		  fprintf(outfp, " max-norm of Scalar ADI-iteration correction ");
		  fprintf(outfp, "for first pde  = (%4f,%4f,%4f)\n", 
			  idmax(1), jdmax(1), kdmax(1));
		  fprintf(outfp, " max-norm of Scalar ADI-iteration correction ");
		  fprintf(outfp, "for second pde = (%4f,%4f,%4f)\n", 
			  idmax(2), jdmax(2), kdmax(2));
		  fprintf(outfp, " max-norm of Scalar ADI-iteration correction ");
		  fprintf(outfp, "for third pde  = (%4f,%4f,%4f)\n", 
			  idmax(3), jdmax(3), kdmax(3));
		  fprintf(outfp, " max-norm of Scalar ADI-iteration correction ");
		  fprintf(outfp, "for fourth pde = (%4f,%4f,%4f)\n", 
			  idmax(4), jdmax(4), kdmax(4));
		  fprintf(outfp, " max-norm of Scalar ADI-iteration correction ");
		  fprintf(outfp, "for fifth pde  = (%4f,%4f,%4f)\n\n", 
			  idmax(5), jdmax(5), kdmax(5));
		  fflush(outfp);
		}
	    }
	  else if ( lnorm == 2 )
	    {
	      l2norm ( isiz1, isiz2, isiz3,
		      nx, ny, nz,
		      &rsd(1,1,1,1), &delunm(1) );

	      if ( ipr == 1 ) 	
		{
		  fprintf(outfp, " RMS-norm of scalar adi-iteration correction ");
		  fprintf(outfp, "for first pde  = %e\n", delunm(1));
		  fprintf(outfp, " RMS-norm of scalar adi-iteration correction ");
		  fprintf(outfp, "for second pde = %e\n", delunm(2));
		  fprintf(outfp, " RMS-norm of scalar adi-iteration correction ");
		  fprintf(outfp, "for third pde  = %e\n", delunm(3));
		  fprintf(outfp, " RMS-norm of scalar adi-iteration correction ");
		  fprintf(outfp, "for fourth pde = %e\n", delunm(4));
		  fprintf(outfp, " RMS-norm of scalar adi-iteration correction ");
		  fprintf(outfp, "for fifth pde  = %e\n\n", delunm(5));
		  fflush(outfp);
		}  
	    }
	}
         
      /*compute the steady-state residuals */
      rhs();

      /*compute the norms of the residuals */
      if ( ( mod ( istep, inorm ) == 0 ) || ( istep == itmax ) ) 
	{
	  if ( lnorm == 1 ) 	
	    {
	      maxnorm ( isiz1, isiz2, isiz3,
		       nx, ny, nz,
		       &imax(1), &jmax(1), &kmax(1),
		       &rsd(1,1,1,1), &rsdnm(1) );
	      
	      if ( ipr == 1 )  
		{
		  fprintf(outfp, "          Initial Residual norms\n\n");
		  
		  fprintf(outfp, " max-norm of steady-state residual ");
		  fprintf(outfp, "for first pde  = (%4f,%4f,%4f)\n", 
			  imax(1), jmax(1), kmax(1));
		  fprintf(outfp, " max-norm of steady-state residual ");
		  fprintf(outfp, "for second pde = (%4f,%4f,%4f)\n", 
			  imax(2), jmax(2), kmax(2));
		  fprintf(outfp, " max-norm of steady-state residual ");
		  fprintf(outfp, "for third pde  = (%4f,%4f,%4f)\n", 
			  imax(3), jmax(3), kmax(3));
		  fprintf(outfp, " max-norm of steady-state residual ");
		  fprintf(outfp, "for fourth pde = (%4f,%4f,%4f)\n", 
			  imax(4), jmax(4), kmax(4));
		  fprintf(outfp, " max-norm of steady-state residual ");
		  fprintf(outfp, "for fifth pde  = (%4f,%4f,%4f)\n", 
			  imax(5), jmax(5), kmax(5));
		  fflush(outfp);
		}
	      }
	  else if ( lnorm == 2 ) 
	    {
	      l2norm ( isiz1, isiz2, isiz3,
		      nx, ny, nz,
		      &rsd(1,1,1,1), &rsdnm(1) );

	      if ( ipr == 1 )  
		{
		  fprintf(outfp, " RMS-norm of steady-state residual for ");
		  fprintf(outfp, "first pde  = %e\n", rsdnm(1));
		  fprintf(outfp, " RMS-norm of steady-state residual for ");
		  fprintf(outfp, "second pde = %e\n", rsdnm(2));
		  fprintf(outfp, " RMS-norm of steady-state residual for ");
		  fprintf(outfp, "third pde  = %e\n", rsdnm(3));
		  fprintf(outfp, " RMS-norm of steady-state residual for ");
		  fprintf(outfp, "fourth pde = %e\n", rsdnm(4));
		  fprintf(outfp, " RMS-norm of steady-state residual for ");
		  fprintf(outfp, "fifth pde  = %e\n\n", rsdnm(5));
		  fflush(outfp);
		}
	    }
	  

	  /*
	  || check the pseudo-time iteration residuals against the 
	  || tolerance levels 
	  */

	  if (( rsdnm(1) < tolrsd(1) ) &&
	      ( rsdnm(2) < tolrsd(2) ) &&
	      ( rsdnm(3) < tolrsd(3) ) &&
	      ( rsdnm(4) < tolrsd(4) ) &&
	      ( rsdnm(5) < tolrsd(5) ) )
	    {
	      fprintf(outfp, "  convergence was achieved after %d pseudo time-steps\n",
		      istep);
	      return;
	    }
	}
      tend = timer( );
      ttotal = tend - tstart;
    }
}
 

/* compute the max-norm of vector v */

maxnorm ( ldmi, ldmj, ldmk, nx, ny, nz, IMAX, JMAX, KMAX, V, VNM )
double V[FOUR_D_SIZE], VNM[5];
int IMAX[5], JMAX[5], KMAX[5];
{
  int m, i, j, k;
  double t1;
  for (m=1;m<=5;m++)
    {
      vnm(m) = -1.0e+10;
    }

  for (k=1; k<=nz;k++) 
    {
      for (j=1; j<=ny;j++) 
	{
	  for (i=1; i<=nx;i++) 
	    {
	      for (m=1;m<=5;m++) 
		{
		  t1 = abs ( v( m, i, j, k ) );
		  if ( vnm(m) < t1 ) 
		    {
		      vnm(m) = t1;
		      imax(m) = i;
		      jmax(m) = j;
		      kmax(m) = k;
		    }
		}
	    }
	}
    }
}


