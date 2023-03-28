#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include "kernels.hpp"
#include "kde_functions.hpp"
#include "../CTvalues_kde.hpp"

/*----------------------------------------------------------------------
 |  Calculations of the caputo derivative
 |  Contains mainly of a scalar and vector version
 |  To Accomadate the dimensions of space, a template the used
 |  The base class contains the precompute for the caputo parameters, which
 |    is the same for all classes
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

namespace kde_functions {


/*----------------------------------------------------------------------
 |  Beta distribution kernels
 -----------------------------------------------------------------------*/
  void kde_gaussian_estimate(const double mean[], const int n_mean, double bandwidth,
    const double x[], const int n_x, double y[])
  {
    int i_start = 0;

    kernels::Gaussian_kernel kernel_gauss;
    kernels::Gamma_kernel    kernel_gamma;

    for (int k = 0; k < n_mean; k++)
    {

      kernel_gauss.set_mean(mean[k], ctv::gauss_bandwitchfactor*bandwidth, n_mean);

      for (int i = i_start; i < n_x; i++)
      {
        if (x[i] > (mean[k] + 2.0*bandwidth))
        {
          break;
        }

        if (x[i] < (mean[k] - 2.0*bandwidth))
        {
          i_start = i_start + 1;
          continue;
        }

      y[i] = y[i] + kernel_gauss.pdf(x[i]);

      }
    }

  }


  void kde_gaussian_bounded_estimate(const double mean[], const int n_mean, double bandwidth,
    const double x[], const int n_x, double y[])
  {
    int i_start = 0;

    kernels::Gaussian_kernel kernel_gauss;
    kernels::Gamma_kernel    kernel_gamma;

    for (int k = 0; k < n_mean; k++)
    {

      if ((mean[k] - ctv::kde_bound_trigger_value * bandwidth) > 0.0)
      {
        kernel_gauss.set_mean(mean[k], ctv::gauss_bandwitchfactor * bandwidth, n_mean);
      }
      else
      {
        kernel_gauss.set_mean(mean[k], ctv::gauss_bandwitchfactor * mean[k], n_mean);
      }

      for (int i = i_start; i < n_x; i++)
      {
        if (x[i] > (mean[k] + 4.0*bandwidth))
        {
          break;
        }

        if (x[i] < (mean[k] - 4.0*bandwidth))
        {
          i_start = i_start + 1;
          continue;
        }

      y[i] = y[i] + kernel_gauss.pdf(x[i]);

      }
    }

  }



}

/*----------------------------------------------------------------------
 |  The precomputed beta and tau values are given here. For now, only
 |  N_p = 9 Prony terms are included (for 500 Fourier series approximations)
 |
 |  More should be imported from the matlab data files at a later time
 -----------------------------------------------------------------------*/