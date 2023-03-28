# File: HModels.pyx
# distutils: language = c++
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
# cython: language_level=3
'''
To be written

'''
import numpy as np
cimport numpy as np
cimport cython

from src.cython.headers.kernel_density_estimation.kde_gaussian cimport (
  kde_gaussian_estimate,
  kde_gaussian_bounded_estimate,
)

# ------------------------------------------------------------------------------
# Residual Functions
#
# COMMENT:
# ------------------------------------------------------------------------------
def KDE_gaussian_approximate(
  double[:] mean,
  double bandwidth,
  double[:] x
):
  n_x    = int(x.shape[0])
  n_mean = int(mean.shape[0])
  cdef np.ndarray[dtype = np.float64_t, ndim=1] y = np.zeros(n_x, dtype=np.float64)
  cdef double[:] c_y = y

  kde_gaussian_estimate(&mean[0], n_mean, bandwidth, &x[0], n_x, &c_y[0])

  return y

def KDE_gaussian_bounded_approximate(
  double[:] mean,
  double bandwidth,
  double[:] x
):
  n_x    = int(x.shape[0])
  n_mean = int(mean.shape[0])
  cdef np.ndarray[dtype = np.float64_t, ndim=1] y = np.zeros(n_x, dtype=np.float64)
  cdef double[:] c_y = y

  kde_gaussian_bounded_estimate(&mean[0], n_mean, bandwidth, &x[0], n_x, &c_y[0])

  return y
