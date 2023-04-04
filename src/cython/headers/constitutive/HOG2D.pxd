# File: HOG2D.pxd
# distutils: language = c++
# cython: language_level=3


""" ----------------------------------------------------------------------------
C++ Source Files
---------------------------------------------------------------------------- """

cimport src.cython.headers.kinematics.tensor_algebra
cimport src.cython.headers.kinematics.kinematics

cdef extern from "src/cpp/constitutive/HOG2D.cpp":
  pass

""" ----------------------------------------------------------------------------
End of Source Files
---------------------------------------------------------------------------- """


# ------------------------------------------------------------------------------
# C++ Header files + exported definitions
# ------------------------------------------------------------------------------

cdef extern from "src/cpp/constitutive/HOG2D.hpp" namespace "constitutive_models":
  cdef cppclass StrucHOG2D:
    StrucHOG2D() except +
    StrucHOG2D(double k1, double k2, double theta, double alpha, double beta, double kip, double kop) except +
    double k1, k2
    double A, B, C
    double m4[4]
    double m6[4]
    double H4[4]
    double H6[4]
    void set_pars(double k1, double k2, double theta, double alpha, double beta, double kip, double kop)
    void stress(double[] args, double[] stress)


  cdef cppclass HOG2D:
    HOG2D() except +
    HOG2D(double k1, double k2, double theta) except +
    double k1, k2
    double m[4]
    void set_pars(double k1, double k2, double theta)
    void stress(double[] args, double[] stress)

