# File: tensor_algebra.pxd
# distutils: language = c++
# cython: language_level=3


""" ----------------------------------------------------------------------------
C++ Source Files
---------------------------------------------------------------------------- """


cdef extern from "src/cpp/kinematics/tensor_algebra.cpp":
  pass


""" ----------------------------------------------------------------------------
End of Source Files
---------------------------------------------------------------------------- """


# ------------------------------------------------------------------------------
# C++ Header files + exported definitions
# ------------------------------------------------------------------------------

cdef extern from "src/cpp/kinematics/tensor_algebra.hpp" namespace "constitutive_models":
  pass

