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

/* common.h */

#ifndef _common_h
#define _common_h

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "input.h"

/*
|| The vanilla code is in Fortran.  We use the following macro transformations 
|| to reuse the original Fortran code with minimal changes. This reduces debugging
|| time.  We rely on the compiler to optimize away the macro transformation 
|| overhead.  We do not to use the f2c output as our baseline because we need to 
|| study the C code, as opposed to just running it.  We feel the f2c output is 
|| pretty ugly and difficult to understand and not suitable for our purpose.
||
||							-Shubu /SSM/
*/

#define TWO_D_SIZE	5 * ISIZ1
#define FOUR_D_SIZE	5*ISIZ1*ISIZ2*ISIZ3
#define THREE_D_SIZE	ISIZ1*ISIZ2*ISIZ3

MAIN int    nx, ny, nz, ii1, ii2, ji1, ji2, ki1, ki2, idum1;
MAIN double dxi, deta, dzeta, tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3;
MAIN double dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, dy5, dz1, dz2, 
       dz3, dz4, dz5, dssp;

#define TWO_D_ARRAY(arrayname, dim1, dim2) \
  (*(arrayname + (dim2 - 1) * 5 + (dim1 - 1)))

/* column major representation */
#define FOUR_D_ARRAY(arrayname, dim1, dim2, dim3, dim4) \
  (*(arrayname + ((dim4-1) * 5 * ISIZ1 * ISIZ2) \
               + ((dim3-1) * 5 * ISIZ1) \
               + ((dim2-1) * 5) \
               + (dim1-1)))

MAIN double U[FOUR_D_SIZE], RSD[FOUR_D_SIZE], FRCT[FOUR_D_SIZE];

#define u(dim1, dim2, dim3, dim4)     FOUR_D_ARRAY(U, dim1, dim2, dim3, dim4)
#define rsd(dim1, dim2, dim3, dim4)   FOUR_D_ARRAY(RSD, dim1, dim2, dim3, dim4)
#define frct(dim1, dim2, dim3, dim4)  FOUR_D_ARRAY(FRCT, dim1, dim2, dim3, dim4)


MAIN int ipr, iout, inorm;
MAIN int itmax, invert;
MAIN double dt, omega, frc, ttotal;

MAIN double TOLRSD[5], RSDNM[5], ERRNM[5];

#define ONE_D_ARRAY(arrayname, dim1) \
  (*(arrayname + (dim1-1)))

#define tolrsd(dim1) 	ONE_D_ARRAY(TOLRSD, dim1)
#define rsdnm(dim1)	ONE_D_ARRAY(RSDNM, dim1)
#define errnm(dim1)	ONE_D_ARRAY(ERRNM, dim1)

#define THREE_D_ARRAY(arrayname, dim1, dim2, dim3) \
  (*(arrayname + ((dim3-1) * ISIZ1 * ISIZ2) \
               + ((dim2-1) * ISIZ1) \
               + (dim1-1)))

MAIN double A[THREE_D_SIZE],
     	    B[THREE_D_SIZE],
            C[THREE_D_SIZE],
            D[THREE_D_SIZE],
            E[THREE_D_SIZE];

#define a(dim1, dim2, dim3) THREE_D_ARRAY(A, dim1, dim2, dim3)
#define b(dim1, dim2, dim3) THREE_D_ARRAY(B, dim1, dim2, dim3)
#define c(dim1, dim2, dim3) THREE_D_ARRAY(C, dim1, dim2, dim3)
#define d(dim1, dim2, dim3) THREE_D_ARRAY(D, dim1, dim2, dim3)
#define e(dim1, dim2, dim3) THREE_D_ARRAY(E, dim1, dim2, dim3)

MAIN double CE[65];	/* ce[5][13] */

#define ce(dim1, dim2) 	(*(CE + ((dim2-1) * 5) + (dim1-1)))

MAIN int setbv(), 
	 setiv(), 
	 erhs(), 
	 adi(), 
	 error(), 
	 maxnorm(), 
	 l2norm(), 
	 exact(), 
	 error(), 
	 jacx(), 
	 jacy(), 
	 jacz(),
	 ninvr(),
	 pintgr(),
	 pinvr(),
	 rhs(),
	 spentax(),	
	 spentax3(),
	 spentay(),
	 spentay3(),
	 spentaz(),
	 spentaz3(),
	 txinvr(),
         tzebar(),
	 verify();
	

MAIN FILE *outfp;

#define mod(mod_d1, mod_d2)	((mod_d1) % (mod_d2))

MAIN int isiz1, isiz2, isiz3;

#define abs(a) (((a)>0) ? (a) : (-1)*(a))

#define max2(a, b) (((a) > (b)) ? (a) : (b))
#define max3(a, b, c) (((a) > (max2((b), (c)))) ? (a) : (max2((b), (c))))
#define max4(a, b, c, d) (((a) > (max3((b), (c), (d)))) ? (a) : (max3((b), (c), (d))))
#define max5(a, b, c, d, e) (((a) > (max4((b), (c), (d), (e)))) ? (a) : (max4((b), (c), (d), (e))))


#define timer()	((double)time(NULL))

#endif /* _common_h */
