# File: residuals.pxd
# distutils: language = c++
# cython: language_level=3


""" ----------------------------------------------------------------------------
C++ Source Files
---------------------------------------------------------------------------- """

cimport src.cython.headers.CTvalues
cimport src.cython.headers.hmodels.simulate

cdef extern from "src/cpp/hmodels/residuals.cpp":
  pass

""" ----------------------------------------------------------------------------
End of Source Files
---------------------------------------------------------------------------- """


# ------------------------------------------------------------------------------
# C++ Header files + exported definitions
# ------------------------------------------------------------------------------

cdef extern from "src/cpp/hmodels/residuals.hpp" namespace "residuals":
  cdef double calculate_hyperE_residual(double[] pars, double[] fiber,
    double[] visco, double Tf,
    double[] args, double[] stress, double[] dt, double[] weights, int[] index, int[] select,
    int n, int dim, int nprot, int skip)

  cdef double calculate_hyperE_residual_scaled(double[] pars, double[] fiber,
    double[] visco, double Tf, double[] Cmax,
    double[] args, double[] stress, double[] dt, double[] weights, double[] deltaCG,
    double[] hysteresis, double[] alphas,
    int[] index, int[] select, int n, int dim, int nprot, int skip)

  cdef double calculate_viscoE_residual_C_M(double[] pars, double[] fiber,
    double[] visco, double Tf,
    double[] args, double[] stress, double[] dt, double[] weights,
    int[] index, int[] select, int n, int dim, int nprot, int skip)

  cdef double calculate_weighted_viscoE_residual_C_M(double[] pars, double[] fiber,
    double[] visco, double Tf,
    double[] args, double[] stress, double[] dt, double[] weights, double[] deltaCG,
    double[] hysteresis, double[] alphas,
    int[] index, int[] select, int n, int dim, int nprot, int skip)

  cdef double calculate_weighted_viscoE_residual_C_M_scaled(double[] pars, double[] fiber,
    double[] visco, double Tf, double[] Cmax,
    double[] args, double[] stress, double[] dt, double[] weights, double[] deltaCG,
    double[] hysteresis, double[] alphas,
    int[] index, int[] select, int n, int dim, int nprot, int skip)



