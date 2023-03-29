# File: caputo.pxd
# distutils: language = c++
# cython: language_level=3


""" ----------------------------------------------------------------------------
C++ Source Files
---------------------------------------------------------------------------- """


cdef extern from "src/cpp/fractional/caputo.cpp":
  pass


""" ----------------------------------------------------------------------------
End of Source Files
---------------------------------------------------------------------------- """


# ------------------------------------------------------------------------------
# C++ Header files + exported definitions
# ------------------------------------------------------------------------------

cdef extern from "src/cpp/fractional/caputo.hpp" namespace "caputo":
  cdef cppclass caputo_init_4:
    caputo_init_4() except +
    caputo_init_4(double alpha, double Tf, double delta) except +
    double alpha, Tf, delta, beta0
    int N
    double betas[9]
    double taus[9]
    double Q[9*4]
    double f_prev[4]
    void set_pars(double alpha, double Tf, double delta)

