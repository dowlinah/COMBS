/* --- appbt.h ---
 *
 * NAS CFD code: APPBT
 *
 * Original FORTRAN version:
 *
 * Author: Sisira Weeratunga
 *         NASA Ames Research Center
 *         (10/25/90)
 *
 * This is the C version of the FORTRAN code developed at NASA Ames Research
 * Center.  We converted the serial code to C as a part of our class project
 * for CS 838-3/ChE 562, which was offered in Spring 1993 by Mark Hill, 
 * Sangtae Kim, and Mary Vernon at the University of Wisconsin-Madison.  
 *
 *					July 20, 1993.
 *
 *					Doug Burger
 *					dburger@mysost.cs.wisc.edu
 *					Computer Sciences Department
 *
 *					Sanjay Mehta
 *					sanjaym@luther.che.wisc.edu
 *					Chemical Engineering Department
 *					
 *					University of Wisconsin-Madison
 *
 */

#define MAXLINE 255

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (double)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (double)min(a,b)
#define dmax(a,b) (double)max(a,b)

#define S1 12 /* Must be >= your problem size */
#define S2 12
#define S3 12

#define XO(src, i) src = GMEM->x_a[i - 1];
#define YO(src, j) src = GMEM->y_a[j - 1];
#define ZO(src, k) src = GMEM->z_a[k - 1];

/* Macros for converting array references to global struct refs */

/*The g's are multiplied by 5 before being passed to save computation */
#define AA(e, f, g) GMEM->a[e + 5 * (f - 1 + g) - 1] 
#define BA(e, f, g) GMEM->b[e + 5 * (f - 1 + g) - 1]
#define CA(e, f, g) GMEM->c[e + 5 * (f - 1 + g) - 1]
#define FRCT GMEM->frct
#define UA GMEM->u
#define RSD GMEM->rsd
#define CE GMEM->ce
#define X_A GMEM->x_a
#define Y_A GMEM->y_a
#define Z_A GMEM->z_a

typedef struct global_struct_ {
	double u[5*S1*S2*S3];
	double rsd[5*S1*S2*S3];
	double frct[5*S1*S2*S3];
	double a[5*5*S1*S2*S3];
	double b[5*5*S1*S2*S3];
	double c[5*5*S1*S2*S3];
	double ce[65];
	int x_a[S1], y_a[S2], z_a[S3];
} Global_Struct;


