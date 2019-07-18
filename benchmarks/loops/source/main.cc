#include "../include/point.h"
#include "../include/integrate.h"

int main()
{
  const unsigned int k=20;
  double u[k][k][k];
  double phi[k][k];
  double omega[k];

  for (unsigned int j=0; j<k; ++j)
    {
      double x = 1.*j/(k-1);
      omega[j] = 1./k;
      phi[0][j] = 1.;
      for (unsigned int i=1; i<k; ++i)
        phi[j][i] = phi[i-1][j]*j/k;
    }

  double value = 0., integral = 0.;

  for (unsigned int i=0; i<1000; ++i)
    {
#ifndef FACTOR
      value += point<k>(u, phi[1]);
      integral += integrate<k,k>(u,phi,omega);
#else
      value += point_fact<k>(u, phi[1]);
      integral += integrate_fact<k,k>(u,phi,omega);
#endif
    }

  std::cout << "Evaluate " << value << ' ' << integral << std::endl;
}

