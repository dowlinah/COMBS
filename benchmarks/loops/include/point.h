#include <iostream>
#include <cmath>


/**
 * Evaluate a function u given as linear combination of shape function
 * values in a single quadrature point
 *
 * \param[in] u \fk^3$ values, coefficients for the linear combination
 * \param[in] phi \f$k$ values of shape functions in a single
 * quadrature point.
 */
template<int k>
double __attribute__ ((noinline)) point(const double u[k][k][k], const double phi[k])
{
  double f = 0.;

  for (int i3 = 0; i3<k; ++i3)
    for (int i2 = 0; i2<k; ++i2)
      for (int i1 = 0; i1<k; ++i1)
        f += u[i1][i2][i3] * phi[i1] * phi[i2] * phi[i3];
  return f;
}

template<int k>
double __attribute__ ((noinline)) point_fact(const double u[k][k][k], const double phi[])
{
  double f3 = 0.;

  for (int i3 = 0; i3<k; ++i3)
    {
      double f2 = 0.;
      for (int i2 = 0; i2<k; ++i2)
        {
          double f1 = 0.;
          for (int i1 = 0; i1<k; ++i1)
            {
              f1 += u[i1][i2][i3] * phi[i1];
            }
          f2 += f1  * phi[i2];
        }
      f3 += f2 * phi[i3];
    }
  return f3;
}

