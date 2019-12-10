#include <stdio.h>
#include <math.h>
#include <stdlib.h>
void heat()
{
	int Nx = 10000;
	int Nt = 1000000;
	double D = 1.14;
	double v = .85;
	double t = 0.1;
	double dt = t/Nt;
	double dx = 1.0/Nx;
	double twopi = 6.283185307179586;

	double *u = malloc(Nx*sizeof(double));
	double *u_next = malloc(Nx*sizeof(double));

	for (int ix = 0; ix<Nx; ix++){
		/* double x = ((double)ix)/Nx; */
		double x = ((double) ix)/(Nx-1);
		u[ix] = sin(twopi*x)-sin(2*twopi*x);
	}
	for(int it = 0;it<Nt;it++){
		for (int ix = 1; ix<Nx-1; ix++){
			/* u_next[ix] = 1.2*u[ix] + 1.4*u[ix-1] + 1.5*u[ix+1]; */
			u_next[ix] = u[ix] + D*dt/(dx*dx)*(u[ix+1]-2*u[ix]+u[ix-1]) - v*dt/dx*(u[ix+1] - u[ix]);
		}
		for (int ix = 0; ix<Nx; ix++){
			u[ix] = u_next[ix];
		}
	}
	printf("%g\n",u[2]);
}
int main()
{
	heat();
	return 0;
}
