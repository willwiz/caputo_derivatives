#!/usr/bin/env python3

import numpy as np
from numpy import ndarray, empty, zeros
from typing import Callable

de = np.cbrt(np.finfo(float).eps)


def hessian(f:Callable, x:ndarray):
  dim = x.size
  f00 = 2.0*f(x)
  hess = empty((dim,dim), dtype=float)
  first_fp = empty((dim, 2), dtype=float)
  h = x.copy()
  h[h<1.0] = 1.0
  h = h*de
  for i in range(dim):
    dh = zeros((dim, ), dtype=float)
    dh[i] = h[i]
    first_fp[i,0] = f(x+dh)
    first_fp[i,1] = f(x-dh)
  for i in range(dim):
    for j in range(dim):
      if i==j:
        hess[i,i] = (first_fp[i,0] - f00 + first_fp[i, 1]) / h[i]**2
      else:
        dh = zeros((dim, ), dtype=float)
        dh[i] = h[i]
        dh[j] = h[j]
        hess[i,j] = f(x+dh) + f(x-dh)
        dh[j] = - h[j]
        hess[i,j] = hess[i,j] - f(x+dh) - f(x-dh)
        hess[i,j] = 0.25 * hess[i,j] / h[i] / h[j]
  return np.around(hess, 6)