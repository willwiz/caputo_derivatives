import os
import numpy as np
from scipy.io import loadmat
from numpy import interp, exp, pi

pack_dir = os.path.dirname(os.path.abspath(__file__))

class make_generic:
  def __init__(self) -> None:
      pass

class caput_init:
  def __init__(self, alpha : float, Tf : float, N : int) -> None:
    freq = 2.0*pi/Tf
    par_file = f"coeffs-opt-refined-100steps-{N}-500.mat"
    carp      = loadmat(os.path.join(pack_dir, "caputo_parm", par_file))
    betam     = carp['betam']
    taum      = carp['taum']
    x         = np.linspace(0.01, 1.0, num=100, endpoint=True)
    self.N     = N
    self.freq  = freq
    self.beta  = np.array([interp(alpha, x, i) * (freq**alpha) for i in betam[:-1]])
    self.beta0 = interp(alpha, x, betam[-1]) * (freq**(alpha - 1.0 ))
    self.tau   = np.array([interp(alpha, x, i) / freq for i in taum])


def caputo_initialize(alpha, Tf, N):
  res       = make_generic()
  freq      = 2*pi/(Tf)
  par_file  = "coeffs-opt-refined-100steps-{:}-500.mat".format(N)
  carp      = loadmat(os.path.join(pack_dir, "caputo_parm", par_file))
  betam     = carp['betam']
  taum      = carp['taum']
  x         = np.linspace(0.01, 1.0, num=100, endpoint=True)
  res.N     = N
  res.freq  = freq
  res.beta  = np.array([interp(alpha, x, i) * (freq**alpha) for i in betam[:-1]])
  res.beta0 = interp(alpha, x, betam[-1]) * (freq**(alpha - 1.0 ))
  res.tau   = np.array([interp(alpha, x, i) / freq for i in taum])
  return res


def caputo_derivative1_array(force: np.ndarray, dt: float, carp: caput_init) -> np.ndarray:
  dim = force.shape
  res = np.zeros(dim)
  store  = np.zeros(carp.N) if (len(dim)==1) else np.zeros((carp.N,dim[1]))
  f_prev = force[0]
  c0  = carp.beta0/dt
  coeff = np.zeros(carp.N)
  for k in range(carp.N):
    coeff[k] = carp.tau[k]/(carp.tau[k] + dt)
  for t in range(dim[0]):
    df = force[t] - f_prev
    f_prev = force[t]
    res[t] = c0 * df 
    for k in range(carp.N):
      store[k] = coeff[k]*(store[k] + carp.beta[k]*df)
      res[t] = res[t] + store[k]
  return res


def caputo_derivative2_array(force: np.ndarray, dt: float, carp: caput_init) -> np.ndarray:
  dim = force.shape
  res = np.zeros(dim)
  store  = np.zeros(carp.N) if (len(dim)==1) else np.zeros((carp.N,dim[1]))
  f_prev = force[0]
  c0  = carp.beta0/dt
  coeffq = np.zeros(carp.N)
  coeffb = np.zeros(carp.N)
  for k in range(carp.N):
    coeffb[k] = exp(-0.5*dt/carp.tau[k])
    coeffq[k] = coeffb[k]*coeffb[k]
  for t in range(dim[0]):
    df = force[t] - f_prev
    f_prev = force[t]
    res[t] = c0 * df 
    for k in range(carp.N):
      store[k] = coeffq[k]*store[k] + carp.beta[k]*coeffb[k]*df 
      res[t] = res[t] + store[k]
  return res


def caputo_diffeq1_array(force: np.ndarray, dt: float, carp: caput_init) -> np.ndarray:
  dim = force.shape
  res = np.zeros(dim)
  store  = np.zeros(carp.N) if (len(dim)==1) else np.zeros((carp.N,dim[1]))
  f_prev = force[0]
  # Precompute
  c0  = carp.beta0/dt
  coeff = np.zeros(carp.N)
  for k in range(carp.N):
    coeff[k] = carp.tau[k]/(carp.tau[k]+dt)
    c0 = c0 + carp.beta[k]*coeff[k]
  c0  = carp.delta*c0
  for t in range(dim[0]):
    res[t] = force[t] + c0 * f_prev
    for k in range(carp.N):
      res[t] = res[t] - carp.delta*coeff[k]*store[k]
    res[t] = res[t]/(1.0 + c0)
    # Update variables
    df = res[t] - f_prev
    f_prev = res[t]
    for k in range(carp.N):
      store[k] = coeff[k]*(store[k] + carp.beta[k]*df)
  return res


    
def caputo_diffeq2_array(force: np.ndarray, dt: float, carp: caput_init) -> np.ndarray:
  dim = force.shape
  res = np.zeros(dim)
  store  = np.zeros(carp.N) if (len(dim)==1) else np.zeros((carp.N,dim[1]))
  f_prev = force[0]
  # Precompute
  c0  = carp.beta0/dt
  coeffb = np.zeros(carp.N)
  coeffq = np.zeros(carp.N)
  for k in range(carp.N):
    coeffb[k] = exp(-0.5*dt/carp.tau[k])
    coeffq[k] = coeffb[k]*coeffb[k]
    c0 = c0 + carp.beta[k]*coeffb[k]
  c0  = carp.delta*c0
  for t in range(dim[0]):
    res[t] = force[t] + c0 * f_prev
    for k in range(carp.N):
      res[t] = res[t] - carp.delta*coeffq[k]*store[k]
    res[t] = res[t]/(1.0 + c0)
    # Update variables
    df = res[t] - f_prev
    f_prev = res[t]
    for k in range(carp.N):
      store[k] = coeffq[k]*store[k] + coeffb[k]*carp.beta[k]*df
  return res