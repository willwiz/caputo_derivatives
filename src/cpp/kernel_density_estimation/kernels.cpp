#define _USE_MATH_DEFINES

#include <math.h>
#include <cmath>
#include <tuple>

#include "kernels.hpp"

/*----------------------------------------------------------------------
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

namespace kernels {


/* ---------------------------------------------------------------------
 |  Beta distribution kernels
  ----------------------------------------------------------------------*/
  Triweight_kernel::Triweight_kernel(/* args */): m_mean{} {}
  Triweight_kernel::Triweight_kernel(double mean, double bandwidth, int n):
    m_mean{mean}, m_1_sigma{1.0/bandwidth}, m_scaling{1.09375/bandwidth/double(n)} {}

  Triweight_kernel::~Triweight_kernel() {}

  void Triweight_kernel::set_mean(double mean, double bandwidth, int n)
  {
    m_mean = mean;
    m_1_sigma = 1.0 / bandwidth;
    m_scaling = 1.09375 * m_1_sigma /double(n);
  }

  double Triweight_kernel::pdf(double x)
  {
    double y = (x - m_mean) * m_1_sigma;

    if (std::abs(y) > 1.0)
    {
      return 0.0;
    }

    double v = (1.0 - y*y);

    return m_scaling * v * v * v;
  }


/* ---------------------------------------------------------------------
 |  Beta distribution kernels
  ----------------------------------------------------------------------*/
  Gaussian_kernel::Gaussian_kernel(/* args */): m_mean{} {}
  Gaussian_kernel::Gaussian_kernel(double mean, double sigma, int n):
    m_mean{mean}, m_1_sigma{1.0/sigma}
  {
    m_scaling = 0.5 * M_SQRT1_2 * M_2_SQRTPI * m_1_sigma / double(n);

  }

  Gaussian_kernel::~Gaussian_kernel() {}

  void Gaussian_kernel::set_mean(double mean, double sigma, int n)
  {
    m_mean = mean;
    m_1_sigma = 1.0 / sigma;
    m_scaling = 0.5 * M_SQRT1_2 * M_2_SQRTPI * m_1_sigma / double(n);
  }

  double Gaussian_kernel::pdf(double x)
  {
    double y = (x - m_mean) * m_1_sigma;

    double v = exp(-0.5*y*y);

    return m_scaling * v;
  }


/* ---------------------------------------------------------------------
 |  Beta distribution kernels
  ----------------------------------------------------------------------*/
  Gamma_kernel::Gamma_kernel(/* args */): m_mean{} {}
  Gamma_kernel::Gamma_kernel(double mean, double sigma, int n):
    m_mean{mean}
  {

    m_b = 10.0 * mean/sigma/sigma;
    m_a = m_b * mean;
    m_scaling = pow(m_b, m_a)/tgamma(m_a);

  }

  Gamma_kernel::~Gamma_kernel() {}

  void Gamma_kernel::set_mean(double mean, double sigma, int n)
  {
    m_mean = mean;

    m_b = 10.0 * mean/sigma/sigma;
    m_a = m_b * mean;
    m_scaling = pow(m_b, m_a)/tgamma(m_a);
  }

  double Gamma_kernel::pdf(double x)
  {
    if (x > 0.0)
    {
      return pow(x, m_a - 1.0) * exp(-m_b*x) * m_scaling;
    }

    return 0.0;
  }

}

