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

/* verify.c */

#define MAIN extern
#include "common.h"

#define xcr(dim1)	XCR[dim1-1]
#define xce(dim1)	XCE[dim1-1]
#define xrr(dim1)	XRR[dim1-1]
#define xre(dim1)	XRE[dim1-1]
  
/* verification routine */
  
verify ( XCR, XCE, xci )
double *XCR, *XCE, xci;
{
  double xri;
  double epsilon;
  int i, j, k, m;
  double XRR[5], XRE[5];
  double tmp;
  
  /* tolerance level */

  epsilon = 1.0e-08;

  if ( ( nx == 12 ) && ( ny == 12 ) && ( nz == 12 ) )
    {
      /* Reference values of RMS-norms of residual, for the (12X12X12) grid, */
      xrr(1) = 2.7470315451339479e-02;
      xrr(2) = 1.0360746705285417e-02;
      xrr(3) = 1.6235745065095532e-02;
      xrr(4) = 1.5840557224455615e-02;
      xrr(5) = 3.4849040609362460e-02;
      
      /* Reference values of RMS-norms of solution error, for the (12X12X12) grid, */
      xre(1) = 2.7289258557377227e-05;
      xre(2) = 1.0364446640837285e-05;
      xre(3) = 1.6154798287166471e-05;
      xre(4) = 1.5750704994480102e-05;
      xre(5) = 3.4177666183390531e-05;
      
      /* Reference value of surface integral, for the (12X12X12) grid, */
      xri = 7.840629330912696200;
      
      /* verification test for residuals */
      for (m=1;m<=5;m++)
	{
	  tmp = fabs ( ( xcr(m) - xrr(m) ) / xrr(m) );
	  if ( tmp > epsilon )
	    {
	      fprintf(outfp, "     VERIFICATION TEST FOR RESIDUALS FAILED\n\n");
	      goto label100;
	    }
	}
      fprintf(outfp, "     VERIFICATION TEST FOR RESIDUALS IS SUCCESSFUL\n\n");
      
      /* verification test for solution error */
      
     label100:     
      for (m=1;m<=5;m++)
	{
	  tmp = fabs ( ( xce(m) - xre(m) ) / xre(m) );
	  if ( tmp > epsilon )
	    {
	      fprintf(outfp, "     VERIFICATION TEST FOR SOLUTION ERRORS FAILED\n\n");
	      goto label200;
	    }
	}
      fprintf(outfp, "     VERIFICATION TEST FOR SOLUTION ERRORS IS SUCCESSFUL\n\n");
      
      /* verification test for surface integral */
      
     label200:
      tmp = fabs ( ( xci - xri ) / xri );
      if ( tmp > epsilon )
	{
	  fprintf(outfp,"     VERIFICATION TEST FOR SURFACE INTEGRAL FAILED\n\n");
	}
      else 
	{
	  fprintf(outfp, "     VERIFICATION TEST FOR SURFACE INTEGRAL IS SUCCESSFUL\n\n\n");
	}

      fprintf(outfp, "          CAUTION\n\n");
      fprintf(outfp, "     REFERENCE VALUES CURRENTLY IN THIS VERIFICATION ROUTINE\n");
      fprintf(outfp, "     ARE VALID ONLY FOR RUNS WITH THE FOLLOWING PARAMETER VALUES:\n\n");
      fprintf(outfp, "     NX = 12;  NY = 12;  NZ = 12\n\n");
      fprintf(outfp, "     ITMAX = 100\n\n");
      fprintf(outfp, "     DT = 1.5e-02\n\n");
      fprintf(outfp, "     CHANGE IN ANY OF THE ABOVE VALUES RENDER THE REFERENCE VALUES\n");
      fprintf(outfp, "     INVALID AND CAUSES A FAILURE OF THE VERIFICATION TEST.\n\n\n");
    }
  else if ( ( nx == 64 ) &&
	   ( ny == 64 ) &&
	   ( nz == 64 ) )
    {
      /* Reference values of RMS-norms of residual, for the (64X64X64) grid, */
      xrr(1) = 2.479982239930019500;
      xrr(2) = 1.127633796436883200;
      xrr(3) = 1.502897788877049100;
      xrr(4) = 1.421781621169517900;
      xrr(5) = 2.129211303513828000;
      
      
      /* Reference values of RMS-norms of solution error, for the (64X64X64) grid, */
      xre(1) = 1.0900140297820550e-04;
      xre(2) = 3.7343951769282091e-05;
      xre(3) = 5.0092785406541633e-05;
      xre(4) = 4.7671093939528255e-05;
      xre(5) = 1.3621613399213001e-04;
      
      /* Reference value of surface integral, for the (64X64X64) grid, */
      xri = 1.208043968503863001e+01;
      
      /* verification test for residuals */
      for (m=1;m<=5;m++)
	{
	  tmp = fabs ( ( xcr(m) - xrr(m) ) / xrr(m) );
	  if ( tmp > epsilon )
	    {
	      fprintf(outfp, "     VERIFICATION TEST FOR RESIDUALS FAILED\n\n");		 
	      goto label400;
	    }
	}
      fprintf(outfp, "     VERIFICATION TEST FOR RESIDUALS IS SUCCESSFUL\n\n");

      /* verification test for solution error */
      
     label400:
      for (m=1;m<=5;m++)
	{
	  tmp = fabs ( ( xce(m) - xre(m) ) / xre(m) );

	  if ( tmp > epsilon )
	    {
	      fprintf(outfp, "     VERIFICATION TEST FOR SOLUTION ERRORS FAILED\n\n");	
	      goto label500;
	    }
	}
      fprintf(outfp, "     VERIFICATION TEST FOR SOLUTION ERRORS IS SUCCESSFUL\n\n");
      
      /* verification test for surface integral */

     label500:
      tmp = fabs ( ( xci - xri ) / xri );
      if ( tmp > epsilon )
	{
	  fprintf(outfp,"     VERIFICATION TEST FOR SURFACE INTEGRAL FAILED\n\n");      
	}
      else
	{
	  fprintf(outfp, "     VERIFICATION TEST FOR SURFACE INTEGRAL IS SUCCESSFUL\n\n");
	}

      fprintf(outfp, "          CAUTION\n\n");
      fprintf(outfp, "     REFERENCE VALUES CURRENTLY IN THIS VERIFICATION ROUTINE\n");
      fprintf(outfp, "     ARE VALID ONLY FOR RUNS WITH THE FOLLOWING PARAMETER VALUES:\n\n");
      fprintf(outfp, "     NX = 64;  NY = 64;  NZ = 64\n\n");
      fprintf(outfp, "     IMAX = 100\n\n");
      fprintf(outfp, "     DT = 1.5e-03\n\n");
      fprintf(outfp, "     CHANGE IN ANY OF THE ABOVE VALUES RENDER THE REFERENCE VALUES\n");
      fprintf(outfp, "     INVALID AND CAUSES A FAILURE OF THE VERIFICATION TEST.\n\n\n");
    }
  else 
    {
      fprintf(outfp, " FOR THE PROBLEM PARAMETERS IN USE NO REFERENCE VALUES ARE PROVIDED\n");
      fprintf(outfp, " IN THE CURRENT VERIFICATION ROUTINE - NO VERIFICATION TEST WAS PERFORMED\n");
    }
}



