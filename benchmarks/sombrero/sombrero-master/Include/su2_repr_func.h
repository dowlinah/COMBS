/***************************************************************************\
* Copyright (c) 2008, Claudio Pica and Agostino Patella                     *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef SUN_REPR_FUNC_H
#define SUN_REPR_FUNC_H

#if (NG!=2) 
#error : Mismatch between NG and representation functions!
#endif
#if (NF!=2) 
#error : Mismatch between NF and representation functions!
#endif



/*
*   _FUND_NORM2 = C(FUND)
*   _FUND_NORM2 = C(REPR)
*/
#define _FUND_NORM2 (+5.000000000000000e-01)
#define _REPR_NORM2 (+5.000000000000000e-01)

/*
*	v COMPLEX matrix 2*2
*	u complex matrix 2*2
*	v is the representation of u, as element of SU(2)
*/
#define _group_represent(v,u) \
	(v).c[0].re = +(u).c[0].re; \
	(v).c[0].im = +(u).c[0].im; \
	(v).c[1].re = +(u).c[1].re; \
	(v).c[1].im = +(u).c[1].im; \
	(v).c[2].re = +(u).c[2].re; \
	(v).c[2].im = +(u).c[2].im; \
	(v).c[3].re = +(u).c[3].re; \
	(v).c[3].im = +(u).c[3].im;


/*
*	m COMPLEX matrix 2*2
*	h real 3-vector
*	m = sum(A) i*h(A)*rT(A)
*/
#define _algebra_represent(m,h) \
	(m).c[0].re = 0.0; \
	(m).c[0].im = +5.000000000000000e-01*(h).c[2]; \
	(m).c[1].re = -5.000000000000000e-01*(h).c[0]; \
	(m).c[1].im = +5.000000000000000e-01*(h).c[1]; \
	(m).c[2].re = +5.000000000000000e-01*(h).c[0]; \
	(m).c[2].im = +5.000000000000000e-01*(h).c[1]; \
	(m).c[3].re = 0.0; \
	(m).c[3].im = -5.000000000000000e-01*(h).c[2];



#ifndef WITH_QUATERNIONS

/*
*	h real 3-vector
*	m COMPLEX matrix 2*2
*	h(A) = -i*ReTr[rT(A).m]/C(r)
*/
#define _algebra_project(h,m) \
	(h).c[0] = -(m).c[1].re+(m).c[2].re; \
	(h).c[1] = +(m).c[1].im+(m).c[2].im; \
	(h).c[2] = -(-(m).c[0].im+(m).c[3].im);


#define _algebra_project_FMAT(h,m) _algebra_project(h,m) 


/*
*	m COMPLEX matrix 2*2
*	h real 3-vector
*	m = sum(A) i*h(A)*T(A)
*/
#define _fund_algebra_represent(m,h) \
	(m).c[0].re = 0.0; \
	(m).c[0].im = +5.000000000000000e-01*(h).c[2]; \
	(m).c[1].re = -5.000000000000000e-01*(h).c[0]; \
	(m).c[1].im = +5.000000000000000e-01*(h).c[1]; \
	(m).c[2].re = +5.000000000000000e-01*(h).c[0]; \
	(m).c[2].im = +5.000000000000000e-01*(h).c[1]; \
	(m).c[3].re = 0.0; \
	(m).c[3].im = -5.000000000000000e-01*(h).c[2];



/*
*	h real 3-vector
*	m COMPLEX matrix 2*2
*	h(A) = -i*ReTr[T(A).m]/C(FUND)
*/
#define _fund_algebra_project(h,m) \
	(h).c[0] = -(m).c[1].re+(m).c[2].re; \
	(h).c[1] = +(m).c[1].im+(m).c[2].im; \
	(h).c[2] = -(-(m).c[0].im+(m).c[3].im);



#else

/* QUATERNION CASE */
/* there is no represented field in this case */

/*
*	h real 3-vector
*	m COMPLEX matrix 2*2
*	h(A) = -i*ReTr[T(A).m]/C(FUND)
*/
#define _fund_algebra_project(h,m) \
	(h).c[0] = 2.*(m).c[2]; \
	(h).c[1] = 2.*(m).c[1]; \
	(h).c[2] = 2.*(m).c[3]

/*
*	h real 3-vector
*	m COMPLEX matrix 2*2
*	h(A) = -i*ReTr[rT(A).m]/C(r)
*/
#define _algebra_project(h,m) _fund_algebra_project(h,m)

/*
*	h real 3-vector
*	m COMPLEX matrix 2*2
*	h(A) = -i*ReTr[rT(A).m]/C(r)
*/
#define _algebra_project_FMAT(h,m) \
	(h).c[0] = -(m).c[1].re+(m).c[2].re; \
	(h).c[1] = +(m).c[1].im+(m).c[2].im; \
	(h).c[2] = -(-(m).c[0].im+(m).c[3].im);



/*
*	m COMPLEX matrix 2*2
*	h real 3-vector
*	m = sum(A) i*h(A)*T(A)
*/
#define _fund_algebra_represent(m,h) \
	(m).c[2] = 0.5*(h).c[0]; \
	(m).c[1] = 0.5*(h).c[1]; \
	(m).c[3] = 0.5*(h).c[2]; \
	(m).c[0] = 0.
    

#endif 

#endif 
