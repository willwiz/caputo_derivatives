# File: kinematics.pxd
# distutils: language = c++
# cython: language_level=3


""" ----------------------------------------------------------------------------
C++ Source Files
---------------------------------------------------------------------------- """


cdef extern from "src/cpp/kinematics/kinematics.cpp":
  pass


""" ----------------------------------------------------------------------------
End of Source Files
---------------------------------------------------------------------------- """


# ------------------------------------------------------------------------------
# C++ Header files + exported definitions
# ------------------------------------------------------------------------------

cdef extern from "src/cpp/kinematics/kinematics.hpp" namespace "kinematics":
  pass

