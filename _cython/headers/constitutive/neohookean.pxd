# File: neohookean.pxd
# distutils: language = c++
# cython: language_level=3


""" ----------------------------------------------------------------------------
C++ Source Files
---------------------------------------------------------------------------- """

# cimport _cython.headers.kinematics.tensor_algebra
# cimport _cython.headers.kinematics.kinematics

cdef extern from "constitutive/neohookean.cpp":
  pass

""" ----------------------------------------------------------------------------
End of Source Files
---------------------------------------------------------------------------- """


# ------------------------------------------------------------------------------
# C++ Header files + exported definitions
# ------------------------------------------------------------------------------

cdef extern from "constitutive/neohookean.hpp" namespace "constitutive_models":
  cdef cppclass neohookean:
    neohookean() except +
    neohookean(double mu) except +
    double mus
    void set_pars(double mu)
    void stress(double[] args, double[] stress)



