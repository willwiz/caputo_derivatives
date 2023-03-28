# File: kde_beta.pxd
# distutils: language = c++
# cython: language_level=3


""" ----------------------------------------------------------------------------
C++ Source Files
---------------------------------------------------------------------------- """

cimport src.cython.headers.kernel_density_estimation.kernels

cdef extern from "src/cpp/kernel_density_estimation/kde_functions.cpp":
  pass


""" ----------------------------------------------------------------------------
End of Source Files
---------------------------------------------------------------------------- """


# ------------------------------------------------------------------------------
# C++ Header files + exported definitions
# ------------------------------------------------------------------------------

cdef extern from "src/cpp/kernel_density_estimation/kde_functions.hpp" namespace "kde_functions":
  void kde_gaussian_estimate(double[] mean, int n_mean, double bandwidth,
    double[] x, int n_x, double[] y)

  void kde_gaussian_bounded_estimate(double[] mean, int n_mean, double bandwidth,
    double[] x, int n_x, double[] y)

