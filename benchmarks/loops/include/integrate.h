#include <iostream>
#include <cmath>

#include <array>
#include <algorithm>

/**
 * Evaluate a function u given as linear combination of shape function
 * values in a single quadrature point
 *
 * \param[in] u \fk^3$ values, coefficients for the linear combination
 * \param[in] phi \f$k$ values of shape functions in \f$k$
 * quadrature points.
 */
template<int nphi, int nquad>
double __attribute__ ((noinline)) integrate(const double u[nphi][nphi][nphi],
                                            const double phi[nquad][nphi],
                                            const double omega[nquad])
{
  double integral = 0.;

  for (int q3 = 0; q3<nquad; ++q3)
    for (int i3 = 0; i3<nphi; ++i3)
      for (int q2 = 0; q2<nquad; ++q2)
        for (int i2 = 0; i2<nphi; ++i2)
          for (int q1 = 0; q1<nquad; ++q1)
            for (int i1 = 0; i1<nphi; ++i1)
              integral += u[i1][i2][i3]
                          * phi[q1][i1] * phi[q2][i2] * phi[q3][i3]
                          * omega[q1] * omega[q2] * omega[q3];
  return integral;
}

template<int nphi, int nquad>
double __attribute__ ((noinline)) integrate_fact(const double u[nphi][nphi][nphi],
                                                 const double phi[nquad][nphi],
                                                 const double omega[nquad])
{
  double integral[nphi];

  for (int q = 0; q<nquad; ++q)
    for (int i = 0; i<nphi; ++i)
      integral[i] += phi[q][i] * omega[q];

  double f3 = 0.;
  for (int i3 = 0; i3<nphi; ++i3)
    {
      double f2 = 0.;
      for (int i2 = 0; i2<nphi; ++i2)
        {
          double f1 = 0.;
          for (int i1 = 0; i1<nphi; ++i1)
            {
              f1 += u[i1][i2][i3] * integral[i1];
            }
          f2 += f1 * integral[i2];
        }
      f3 += f2 * integral[i3];
    }
  return f3;
}

