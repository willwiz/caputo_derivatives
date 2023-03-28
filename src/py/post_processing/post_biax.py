#!/usr/bin/env

from typing import Tuple
from numpy import ndarray, array, mean, reshape, sqrt, zeros, einsum
from numpy.linalg import norm
from scipy.linalg import lstsq


dfdr = 0.25*array([-1,-1,1,1],dtype=float)
dfds = 0.25*array([-1,1,-1,1],dtype=float)
gradv = 0.25*array([[-1,-1,1,1],[-1,1,-1,1]],dtype=float)

ex = array([1,0],dtype=float)
ey = array([0,1],dtype=float)
A  = zeros((4,3), dtype=float)


def dudr(x:ndarray) -> Tuple[float,float]:
  return x@dfdr, x@dfds


def grad_tensor(x:ndarray)->ndarray:
  return einsum('ij,kj->ik',x,gradv)


class deformation_gradient:
  def __init__(self, x0:ndarray, y0:ndarray) -> None:
    dXdr = x0 @ dfdr
    dXds = x0 @ dfds
    dYdr = y0 @ dfdr
    dYds = y0 @ dfds
    det  = dXdr*dYds - dXds*dYdr
    self.ref_tensor = array([[dYds, -dXds],[-dYdr,dXdr]],dtype=float)/det
  def deformation_gradient(self, coord:ndarray)->ndarray:
    grad = einsum('mij,kj->mik', coord, gradv)
    DG   = einsum('mij,jk->mik', grad, self.ref_tensor)
    return DG


def find_normals(x,y) -> ndarray:
  n1 = array([y[3]-y[2],x[2]-x[3]])
  n1 = n1/norm(n1)
  n2 = array([y[1]-y[3],x[3]-x[1]])
  n2 = n2/norm(n2)
  n3 = array([y[1]-y[0],x[0]-x[1]])
  n3 = n3/norm(n3)
  n4 = array([y[1]-y[2],x[2]-x[1]])
  n4 = n4/norm(n4)
  return n1,n2,n3,n4


def regression_tensor(n1:ndarray, n2:ndarray, n3:ndarray, n4:ndarray) -> ndarray:
  return 0.5*array([
    [0,0,0,0,0,0,n1[0],n1[1],0,n1[0],n1[1],0],
    [0,0,0,0,n2[0],n2[1],0,0,0,0,n2[0],n2[1]],
    [n3[0],n3[1],0,n3[0],n3[1],0,0,0,0,0,0,0],
    [0,n4[0],n4[1],0,0,0,0,n4[0],n4[1],0,0,0]
  ])


def stress_regression(x,y,f1,f2)-> ndarray:
  n1,n2,n3,n4 = find_normals(x,y)
  A = regression_tensor(n1,n2,n3,n4)
  b = array([f1,f2,f1,f2])
  sigmas = lstsq(A,b,lapack_driver='gelsy', check_finite=False)[0]
  return mean(reshape(sigmas, (-1,3)), axis=0)


def stress_homogenous(tFinv, f1, f2)->ndarray:
  n1 = ex @ tFinv
  n2 = ey @ tFinv
  A[0,  :2] = n1
  A[1, 1: ] = n1
  A[2,  :2] = n2
  A[3, 1: ] = n2
  b = array([f1, 0, 0, f2], dtype=float)
  return lstsq(A,b,lapack_driver='gelsy', check_finite=False)[0]


def mat_vec_contraction(tensor_list, vec):
  mi = einsum("mij,j->mi", tensor_list, vec)
  m  = einsum("mi,mi->m",  mi, mi)
  return sqrt(m)
