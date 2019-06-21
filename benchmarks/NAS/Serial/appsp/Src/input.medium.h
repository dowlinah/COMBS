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

/* input.medium.h */

#ifndef _input_h
#define _input_h

#define ISIZ1 		64
#define ISIZ2		64
#define ISIZ3		64
#define c1		1.4
#define c2		0.4
#define c3		0.1
#define c4		1.0
#define c5		1.4

#define IPR		1
#define INORM		1
#define ITMAX		400
#define DT		0.0015
#define OMEGA		1.2
#define	Tolrsd		1.0e-12
#define NX		64
#define NY		64
#define NZ		64

#endif /* _input_h */

