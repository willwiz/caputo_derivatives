# File: kernels.pxd
# distutils: language = c++
# cython: language_level=3


""" ----------------------------------------------------------------------------
C++ Source Files
---------------------------------------------------------------------------- """


cdef extern from "src/cpp/kernel_density_estimation/kernels.cpp":
  pass


""" ----------------------------------------------------------------------------
End of Source Files
---------------------------------------------------------------------------- """


# ------------------------------------------------------------------------------
# C++ Header files + exported definitions
# ------------------------------------------------------------------------------

cdef extern from "src/cpp/kernel_density_estimation/kernels.hpp" namespace "kde_kernels":
  pass

