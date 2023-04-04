#!/usr/bin/env python3

from .mittag_leffler.mittag_leffler import ml
from scipy.special import gamma

def force_1(a,t):
  return 1.0

def force_t(a,t):
  return t / 4

def force_p(a,t):
  n = len(a)
  tk = 1.0
  res = 0.0
  for i in range(n):
    tk = tk*t
    res = res + a[i]*tk
  return res

def frac_t_anal(par, a, d, t):
  return t**(1.0-a) / gamma(2.0-a) /4

def frac_tk_anal(par, a, k, t):
  return t**(k-a) * gamma(k+1)/ gamma(k+1-a)

def diff_1_anal(par, a, d, t):
  return t**a / d * ml(-t**a / d, a, a+1.0)

def diff_t_anal(par, a, d, t):
  return t**(a+1.0) / d * ml(-t**a / d, a, a+2.0) / 4

def diff_frac_t_anal(par, a, d, t):
  return t / d * ml(-t**a / d, a, 2.0) / 4

def diff_frac_p_anal(par, a, d, t):
  tk = 1.0
  res = 0.0
  for i in range(len(par)):
    tk = tk * t
    res = res + par[i] * tk * gamma(2.0 + i) / d * ml(-t**a /d, a, i + 2.0)
  return res

def diff_frac_p1_anal(par, a, d, t):
  tk = 1.0
  res = 0.0
  for i in range(len(par)):
    tk = tk * t
    res = res + par[i] * tk * gamma(2.0 + i) / d * ml(-t**a /d, a, i + 2.0)
  return res

def diff_frac_p2_anal(par, a, d, t):
  tk = 1.0
  res = 0.0
  for i in range(len(par)):
    tk = tk * t
    res = res + par[i] * tk * gamma(2.0 + i) / d * ml(-t**a /d, a, i + 2.0)
  return res

def diff_frac_p3_anal(par, a, d, t):
  tk = 1.0
  res = 0.0
  for i in range(len(par)):
    tk = tk * t
    res = res + par[i] * tk * gamma(2.0 + i) / d * ml(-t**a /d, a, i + 2.0)
  return res
