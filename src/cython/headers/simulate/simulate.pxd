# File: simulate.pxd
# distutils: language = c++
# cython: language_level=3


""" ----------------------------------------------------------------------------
C++ Source Files
---------------------------------------------------------------------------- """
cimport src.cython.headers.CompileTimeValues
cimport src.cython.headers.constitutive.neohookean
cimport src.cython.headers.constitutive.HOG2D
cimport src.cython.headers.fractional.caputo

cdef extern from "src/cpp/simulate/simulate.cpp":
  pass

""" ----------------------------------------------------------------------------
End of Source Files
---------------------------------------------------------------------------- """


# ------------------------------------------------------------------------------
# C++ Header files + exported definitions
# ------------------------------------------------------------------------------

cdef extern from "src/cpp/simulate/simulate.hpp" namespace "sim":

  cdef void get_model_parameters(double[] pars, double[] fiber, double[] visco, double Tf,
    double[11] scl_pars)

  cdef void get_model_parameters_scaled(double[] pars, double[] fiber, double[] visco, double Tf,
    double[] Cmax, double[11] scl_pars)

  cdef void he_simulate(double[] pars, double[] fiber, double[] caputo, double Tf,
    double[] args, double[] dt, double[] stress, int n)

  cdef void he_simulate_scaled(double[] pars, double[] fiber, double[] caputo, double Tf, double[] Cmax,
    double[] args, double[] dt, double[] stress, int n)

  cdef void caputo_simulate_C_M(double[] pars, double[] fiber, double[] caputo, double Tf,
    double[] args, double[] dt, double[] stress, int n)

  cdef void caputo_simulate_C_M_scaled(double[] pars, double[] fiber, double[] caputo,
    double Tf, double[] Cmax, double[] args, double[] dt, double[] stress, int n)

  cdef void caputo_simulate_all(double[] pars, double[] fiber, double[] caputo, double Tf,
    double[] args, double[] dt, double[] stress, int n)



