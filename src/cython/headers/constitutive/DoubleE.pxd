# File: DoubleE.pxd
# distutils: language = c++
# cython: language_level=3


""" ----------------------------------------------------------------------------
C++ Source Files
---------------------------------------------------------------------------- """

cimport src.cython.headers.kinematics.tensor_algebra
cimport src.cython.headers.kinematics.kinematics

cdef extern from "src/cpp/constitutive/DoubleE.cpp":
  pass


""" ----------------------------------------------------------------------------
End of Source Files
---------------------------------------------------------------------------- """


# ------------------------------------------------------------------------------
# C++ Header files + exported definitions
# ------------------------------------------------------------------------------

cdef extern from "src/cpp/constitutive/DoubleE.hpp" namespace "constitutive_models":
  pass

