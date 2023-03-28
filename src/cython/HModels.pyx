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

from src.cython.headers.constitutive.caputo cimport caputo_init_4
from src.cython.headers.constitutive.neohookean cimport neohookean
from src.cython.headers.constitutive.HOG2D cimport (
  StrucHOG2D,
  HOGDouble2D,
  HOG2D,
)
from src.cython.headers.hmodels.simulate cimport (
  get_model_parameters,
  get_model_parameters_scaled,
  he_simulate,
  he_simulate_scaled,
  caputo_simulate_C_M,
  caputo_simulate_C_M_scaled,
)
from src.cython.headers.hmodels.residuals cimport (
  calculate_hyperE_residual,
  calculate_hyperE_residual_scaled,
  calculate_viscoE_residual_C_M,
  calculate_weighted_viscoE_residual_C_M,
  calculate_weighted_viscoE_residual_C_M_scaled,
)


# ------------------------------------------------------------------------------
# Caputo Derivative
#
# COMMENT:
# ------------------------------------------------------------------------------

cdef class caputo_initialize:
  cdef caputo_init_4 carp
  def __init__(self, double alpha, double Tf, double delta):
    self.carp = caputo_init_4(alpha, Tf, delta)
  def set_pars(self,  double alpha, double Tf, double delta):
    self.carp.set_pars(alpha, Tf, delta)
  @property
  def beta0(self):
    return self.carp.beta0
  @property
  def betas(self):
    return self.carp.betas
  @property
  def taus(self):
    return self.carp.taus


# ------------------------------------------------------------------------------
# Neohookean
#
# COMMENT:
# ------------------------------------------------------------------------------

cdef class NeoHookean2D:
  cdef neohookean model
  def __init__(self, double mu):
    self.model = neohookean(mu)
  def set_pars(self, double mu):
    self.model.set_pars(mu)
  def stress(self, np.ndarray[np.float64_t, ndim=1] args):
    cdef np.ndarray[dtype = np.float64_t, ndim=1] sigma = np.zeros(4, dtype=np.float64)
    self.model.stress(&args[0], &sigma[0])
    return sigma


# ------------------------------------------------------------------------------
# Fiber Models
#
# COMMENT:
# ------------------------------------------------------------------------------

cdef class HOGstruc2D:
  cdef StrucHOG2D model
  def __init__(self, double k1, double k2, double theta, double alpha, double beta, double kip, double kop):
    self.model = StrucHOG2D(k1, k2, theta, alpha, beta, kip, kop)
  def set_pars(self, double k1, double k2, double theta, double alpha, double beta, double kip, double kop):
    self.model.set_pars(k1, k2, theta, alpha, beta, kip, kop)
  def stress(self, np.ndarray[np.float64_t, ndim=1] args):
    cdef np.ndarray[dtype = np.float64_t, ndim=1] sigma = np.zeros(4, dtype=np.float64)
    self.model.stress(&args[0], &sigma[0])
    return sigma


cdef class pyHOGDouble2D:
  cdef HOGDouble2D model
  def __init__(self, double k1, double k2, double theta, double alpha):
    self.model = HOGDouble2D(k1, k2, theta, alpha)
  def set_pars(self, double k1, double k2, double theta, double alpha):
    self.model.set_pars(k1, k2, theta, alpha)
  def stress(self, np.ndarray[np.float64_t, ndim=1] args):
    cdef np.ndarray[dtype = np.float64_t, ndim=1] sigma = np.zeros(4, dtype=np.float64)
    self.model.stress(&args[0], &sigma[0])
    return sigma


cdef class pyHOG2D:
  cdef HOG2D model
  def __init__(self, double k1, double k2, double theta):
    self.model = HOG2D(k1, k2, theta)
  def set_pars(self, double k1, double k2, double theta):
    self.model.set_pars(k1, k2, theta)
  def stress(self, np.ndarray[np.float64_t, ndim=1] args):
    cdef np.ndarray[dtype = np.float64_t, ndim=1] sigma = np.zeros(4, dtype=np.float64)
    self.model.stress(&args[0], &sigma[0])
    return sigma


# ------------------------------------------------------------------------------
# Simulate
#
# COMMENT:
# ------------------------------------------------------------------------------

def GetScaledParameters_cpp(double[:] pars, double[:] fiber, double[:] visco, double Tf):
  cdef np.ndarray[dtype = np.float64_t, ndim=1] scaled_pars = np.zeros((11), dtype=np.float64)
  cdef double[:] parsview = scaled_pars
  get_model_parameters(&pars[0], &fiber[0], &visco[0], Tf, &parsview[0])
  return scaled_pars

def GetScaledParametersScaled_cpp(double[:] pars, double[:] fiber, double[:] visco, double Tf,
    double[:] Cmax):
  cdef np.ndarray[dtype = np.float64_t, ndim=1] scaled_pars = np.zeros((11), dtype=np.float64)
  cdef double[:] parsview = scaled_pars
  get_model_parameters_scaled(&pars[0], &fiber[0], &visco[0], Tf, &Cmax[0], &parsview[0])
  return scaled_pars

def HESimulation(
  double[:] pars,
  double[:] fiber,
  double[:] caputo,
  double Tf,
  double[:,:] args,
  double[:] dt
):
  n = int(dt.shape[0])
  cdef np.ndarray[dtype = np.float64_t, ndim=2] sigma = np.zeros((n, 4), dtype=np.float64)
  cdef double[:, :] stress = sigma
  he_simulate(&pars[0], &fiber[0], &caputo[0], Tf, &args[0][0], &dt[0], &stress[0][0], n)
  return sigma


def HESimulation_Scaled(double[:] pars, double[:] fiber, double[:] caputo,
    double Tf, double[:] Cmax, double[:,:] args, double[:] dt):
  # determine the number of time steps
  n = int(dt.shape[0])
  # create a numpy array to store the results, get memoryview to pass by reference
  # need to make sure that the numpy array is C contiguous
  cdef np.ndarray[dtype = np.float64_t, ndim=2] sigma = np.zeros((n, 4), dtype=np.float64)
  cdef double[:, :] stress = sigma
  he_simulate_scaled(&pars[0], &fiber[0], &caputo[0], Tf, &Cmax[0], &args[0][0], &dt[0], &stress[0][0], n)
  return sigma


def CaputoSimulation_C_M(
  double[:] pars,
  double[:] fiber,
  double[:] caputo,
  double Tf,
  double[:,:] args,
  double[:] dt,
):
  n = int(dt.shape[0])
  cdef np.ndarray[dtype = np.float64_t, ndim=2] sigma = np.zeros((n, 4), dtype=np.float64)
  cdef double[:, :] stress = sigma
  caputo_simulate_C_M(&pars[0], &fiber[0], &caputo[0], Tf, &args[0][0], &dt[0], &stress[0][0], n)
  return sigma

def CaputoSimulation_C_M_Scaled(
  double[:] pars,
  double[:] fiber,
  double[:] caputo,
  double Tf,
  double[:] Cmax,
  double[:,:] args,
  double[:] dt,
):
  n = int(dt.shape[0])
  cdef np.ndarray[dtype = np.float64_t, ndim=2] sigma = np.zeros((n, 4), dtype=np.float64)
  cdef double[:, :] stress = sigma
  caputo_simulate_C_M_scaled(&pars[0], &fiber[0], &caputo[0], Tf, &Cmax[0], &args[0][0], &dt[0], &stress[0][0], n)
  return sigma


# ------------------------------------------------------------------------------
# Residual Functions
#
# COMMENT:
# ------------------------------------------------------------------------------

def ResidualFun_HE(
  double[:] pars,
  double[:] fiber,
  double[:] caputo,
  double Tf,
  double[:,:] strain,
  double[:,:] stress,
  double[:] dt,
  double[:,:] weights,
  int[:] index,
  int[:] select,
  int skip = 1
):
  n     = int(strain.shape[0])
  dim   = int(strain.shape[1])
  nprot = int(select.shape[0])
  return calculate_hyperE_residual(&pars[0], &fiber[0], &caputo[0], Tf,
    &strain[0][0], &stress[0][0], &dt[0], &weights[0][0], &index[0], &select[0],
    n, dim, nprot, int(skip))


def ResidualFun_HE_Scaled(
  double[:] pars,
  double[:] fiber,
  double[:] visco,
  double Tf,
  double[:] Cmax,
  double[:,:] strain,
  double[:,:] stress,
  double[:] dt,
  double[:,:] weights,
  int[:] index,
  int[:] select,
  int skip = 1
):
  n     = int(strain.shape[0])
  dim   = int(strain.shape[1])
  nprot = int(select.shape[0])
  return calculate_hyperE_residual(&pars[0], &fiber[0], &visco[0], Tf,
    &strain[0][0], &stress[0][0], &dt[0], &weights[0][0], &index[0], &select[0],
    n, dim, nprot, int(skip))


def ResidualFun_VE_C_M(
  double[:] pars,
  double[:] fiber,
  double[:] caputo,
  double Tf,
  double[:,:] strain,
  double[:,:] stress,
  double[:] dt,
  double[:,:] weights,
  int[:] index,
  int[:] select,
  int skip = 1
):
  n     = int(strain.shape[0])
  dim   = int(strain.shape[1])
  nprot = int(select.shape[0])
  return calculate_viscoE_residual_C_M(&pars[0], &fiber[0], &caputo[0], Tf,
    &strain[0][0], &stress[0][0], &dt[0], &weights[0][0], &index[0], &select[0],
    n, dim, nprot, int(skip))


def ResidualFun_Weighted_VE_C_M(
  double[:] pars,
  double[:] fiber,
  double[:] caputo,
  double Tf,
  double[:,:] strain,
  double[:,:] stress,
  double[:] dt,
  double[:,:] weights,
  double[:,:] deltaCG,
  double[:,:] hysteresis,
  double[:] alpha,
  int[:] index,
  int[:] select,
  int skip = 1
):
  n     = int(strain.shape[0])
  dim   = int(strain.shape[1])
  nprot = int(select.shape[0])
  return calculate_weighted_viscoE_residual_C_M(&pars[0], &fiber[0], &caputo[0], Tf,
    &strain[0][0], &stress[0][0], &dt[0], &weights[0][0],
    &deltaCG[0][0], &hysteresis[0][0], &alpha[0],
    &index[0], &select[0], n, dim, nprot, int(skip))


def ResidualFun_Weighted_VE_C_M_Scaled(
  double[:] pars,
  double[:] fiber,
  double[:] caputo,
  double Tf,
  double[:] Cmax,
  double[:,:] strain,
  double[:,:] stress,
  double[:] dt,
  double[:,:] weights,
  double[:,:] deltaCG,
  double[:,:] hysteresis,
  double[:] alpha,
  int[:] index,
  int[:] select,
  int skip = 1
):
  n     = int(strain.shape[0])
  dim   = int(strain.shape[1])
  nprot = int(select.shape[0])
  return calculate_weighted_viscoE_residual_C_M_scaled(&pars[0], &fiber[0], &caputo[0], Tf, &Cmax[0],
    &strain[0][0], &stress[0][0], &dt[0], &weights[0][0],
    &deltaCG[0][0], &hysteresis[0][0], &alpha[0],
    &index[0], &select[0], n, dim, nprot, int(skip))