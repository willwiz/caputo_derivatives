import os
import numpy as np
import opt_einsum.contract as einsum
from typing import Tuple, Optional
from scipy.io import loadmat
from numpy import interp, exp, pi, ndarray, zeros

pack_dir = os.path.dirname(os.path.abspath(__file__))


class caputo_init():
  def __init__(self, alpha:float, Tf:float, N:int, delta:float=0.0, dim:Optional[int]=None) -> None:
    freq      = 2*pi/(Tf)
    par_file  = "coeffs-opt-refined-100steps-{:}-500.mat".format(N)
    carp      = loadmat(os.path.join(pack_dir, "caputo_parm", par_file))
    betam     = carp['betam']
    taum      = carp['taum']
    x         = np.linspace(0.01, 1.0, num=100, endpoint=True)
    self.alpha = alpha
    self.N     = N
    self.Tf    = Tf
    self.delta = delta
    self.beta0 = interp(alpha, x, betam[-1]) * (freq**(alpha - 1.0 ))
    self.beta  = np.array([interp(alpha, x, i) * (freq**alpha) for i in betam[:-1]])
    self.tau   = np.array([interp(alpha, x, i) / freq for i in taum])
    if dim is None:
      self.Q = zeros(N)
      self.f_prev = 0.0
    else:
      self.Q = zeros((N,dim))
      self.f_prev = zeros(dim)



def caputo_derivative1_iter(fn : ndarray, dt : float, carp:caputo_init) -> Tuple[ndarray,caputo_init]:
  df   = (fn - carp.f_prev)
  ek = carp.tau / (carp.tau + dt)
  for k in range(carp.N):
    carp.Q[k,:] = ek[k] * (carp.Q[k,:] + carp.beta[k] * df)
  v = (carp.beta0 / dt) * df + einsum('ki->i', carp.Q)
  carp.f_prev = fn
  return v, carp


def caputo_derivative2_iter(fn : ndarray, dt : float, carp:caputo_init) -> Tuple[ndarray,caputo_init]:
  df   = (fn - carp.f_prev)
  ek   = exp( - 0.5 * dt / carp.tau)
  e2   = ek * ek
  for k in range(carp.N):
    carp.Q[k,:] = e2[k] * carp.Q[k,:] +  carp.beta[k] * ek[k] * df
  v = (carp.beta0 / dt) * df + einsum('ki->i', carp.Q)
  # v = (carp.beta0 / dt) * df
  carp.f_prev = fn
  return v, carp


def diffeq_approx1_iter(fn : ndarray, dt : float, carp:caputo_init) -> Tuple[ndarray,caputo_init]:
  K0 = carp.beta0 / dt
  ek = carp.tau / (carp.tau + dt)
  K0 = carp.delta * (K0 + einsum('k,k->', carp.beta, ek))
  v  = fn - carp.delta * einsum('k,ki->k', ek, carp.Q)
  v  = (v + K0 * carp.f_prev)/ (1.0 + K0)
  # Updates
  df = v - carp.f_prev
  carp.f_prev = v
  for k in range(carp.N):
    carp.Q[k,:] = ek[k] * (carp.Q[k,:] + carp.beta[k]*df)
  return v, carp


def diffeq_approx2_iter(fn : ndarray, dt : float, carp:caputo_init) -> Tuple[ndarray,caputo_init]:
  K0 = carp.beta0 / dt
  ek = exp( - 0.5 * dt / carp.tau)
  e2 = ek * ek
  K0 = carp.delta * (K0 + einsum('k,k->', carp.beta, ek))
  v  = fn - carp.delta * einsum('k,ki->k', e2, carp.Q)
  v  = (v + K0 * carp.f_prev)/ (1.0 + K0)
  # Updates
  df = v - carp.f_prev
  carp.f_prev = v
  for k in range(carp.N):
    carp.Q[k,:] = e2[k] * carp.Q[k,:] + carp.beta[k] * ek[k] * df
  return v, carp


