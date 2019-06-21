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

/* main.c */

#define MAIN 
#include "common.h"

/*
|| driver for the performance evaluation of the solver for
|| five coupled, nonlinear partial differential equations.
*/

main()
{
  int i;

  /* Initialize variables */
  ipr   = IPR;
  inorm = INORM;
  itmax = ITMAX;
  dt 	= DT;
  omega = OMEGA;

  for (i=1;i<=5;i++)
    {
      tolrsd(i) = Tolrsd;
    }

  nx = NX;
  ny = NY;
  nz = NZ;
  
  isiz1 = ISIZ1;
  isiz2 = ISIZ2;
  isiz3 = ISIZ3;

  outfp = fopen("output.data", "w");

  if (nx < 5 || ny < 5 || nz < 5) 
    {
      fprintf(outfp, "Problem size is too small");
      fprintf(outfp, "- set each of NX, NY and NZ at least equal to 5");
    }

  if (nx > ISIZ1 || ny > ISIZ2 || nz > ISIZ3)
    {
      fprintf(outfp, "Problem size is too large - NX, NY and NZ should");
      fprintf(outfp, "be less than or equal to ISIZ1, ISIZ2 and ISIZ3");
      fprintf(outfp, "respectively\n");
    }

  dxi = 1.0/ (nx - 1);
  deta = 1.0/ (ny - 1);
  dzeta = 1.0/ (nz - 1);
  
  tx1 = 1.0/ (dxi * dxi);
  tx2 = 1.0/ (dxi * 2);
  tx3 = 1.0/ dxi;
  
  ty1 = 1.0/ (deta * deta);
  ty2 = 1.0/ (deta * 2);
  ty3 = 1.0/ deta;

  tz1 = 1.0/ (dzeta * dzeta);
  tz2 = 1.0/ (dzeta * 2);
  tz3 = 1.0/ dzeta;
  
  ii1 = 2;
  ii2 = nx - 1;
  ji1 = 2;
  ji2 = ny - 2;
  ki1 = 3;
  ki2 = nz - 1;

  /* diffusion coefficients */

  dx1 = .75;
  dx2 = dx1;
  dx3 = dx1;
  dx4 = dx1;
  dx5 = dx1;
  
  dy1 = .75;
  dy2 = dy1;
  dy3 = dy1;
  dy4 = dy1;
  dy5 = dy1;
  
  dz1 = 1.;
  dz2 = dz1;
  dz3 = dz1;
  dz4 = dz1;
  dz5 = dz1;
  
  /* fourth difference dissipation */

  dssp  = max3(dx1, dy1, dz1)/4;

  /* coefficients of the exact solution to the first pde */

  ce(1,1) = 2.0;
  ce(1,2) = 0.0;
  ce(1,3) = 0.0;
  ce(1,4) = 4.0;
  ce(1,5) = 5.0;
  ce(1,6) = 3.0;
  ce(1,7) = 5.0e-01;
  ce(1,8) = 2.0e-02;
  ce(1,9) = 1.0e-02;
  ce(1,10) = 3.0e-02;
  ce(1,11) = 5.0e-01;
  ce(1,12) = 4.0e-01;
  ce(1,13) = 3.0e-01;

  /* coefficients of the exact solution to the second pde */

  ce(2,1)  = 1.0;
  ce(2,2)  = 0.0;
  ce(2,3)  = 0.0;
  ce(2,4)  = 0.0;
  ce(2,5)  = 1.0;
  ce(2,6)  = 2.0;
  ce(2,7)  = 3.0;
  ce(2,8)  = 1.0e-02;
  ce(2,9)  = 3.0e-02;
  ce(2,10) = 2.0e-02;
  ce(2,11) = 4.0e-01;
  ce(2,12) = 3.0e-01;
  ce(2,13) = 5.0e-01;
  
  /* coefficients of the exact solution to the third pde */

  ce(3,1)  = 2.0;
  ce(3,2)  = 2.0;
  ce(3,3)  = 0.0;
  ce(3,4)  = 0.0;
  ce(3,5)  = 0.0;
  ce(3,6)  = 2.0;
  ce(3,7)  = 3.0;
  ce(3,8)  = 4.0e-02;
  ce(3,9)  = 3.0e-02;
  ce(3,10) = 5.0e-02;
  ce(3,11) = 3.0e-01;
  ce(3,12) = 5.0e-01;
  ce(3,13) = 4.0e-01;
  
  /* coefficients of the exact solution to the fourth pde */

  ce(4,1)  = 2.0;
  ce(4,2)  = 2.0;
  ce(4,3)  = 0.0;
  ce(4,4)  = 0.0;
  ce(4,5)  = 0.0;
  ce(4,6)  = 2.0;
  ce(4,7)  = 3.0;
  ce(4,8)  = 3.0e-02;
  ce(4,9)  = 5.0e-02;
  ce(4,10) = 4.0e-02;
  ce(4,11) = 2.0e-01;
  ce(4,12) = 1.0e-01;
  ce(4,13) = 3.0e-01;
  
  /* coefficients of the exact solution to the fifth pde */

  ce(5,1)  = 5.0;
  ce(5,2)  = 4.0;
  ce(5,3)  = 3.0;
  ce(5,4)  = 2.0;
  ce(5,5)  = 1.0e-01;
  ce(5,6)  = 4.0e-01;
  ce(5,7)  = 3.0e-01;
  ce(5,8)  = 5.0e-02;
  ce(5,9)  = 4.0e-02;
  ce(5,10) = 3.0e-02;
  ce(5,11) = 1.0e-01;
  ce(5,12) = 3.0e-01;
  ce(5,13) = 2.0e-01;


  /* set the boundary values for dependent variables */

  setbv(); 	

  /* set the initial values for dependent variables */

  setiv(); 	

  /* compute the forcing term based on prescribed exact solution */

  erhs();  	

  /* perform scalar approximate factorization iterations */

  adi();   	

  /* compute the solution error */

  error(); 	

 /* compute the surface integral */

  pintgr();	

  /* verification test */

  verify(&rsdnm(1), &errnm(1), frc); 	

  /* print the CPU time */

  fprintf(outfp, "     Total CPU time = %e Sec.\n", ttotal); 

  fclose(outfp);
} /* main */

