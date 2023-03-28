#!/usr/bin/env python3
import os, sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from numpy import empty, ascontiguousarray, mean, empty_like, sqrt, zeros_like, absolute
from scipy import interpolate, stats, linalg
try:
  from src.py.post_processing.hermite_poly import hermite_fit, hermite_interpolation, get_hermite_matrix
  from src.py.post_processing.post_biax import *
except:
  from .hermite_poly import hermite_fit, hermite_interpolation, get_hermite_matrix
  from .post_biax import *
import opt_einsum as oe
from time import perf_counter
# import opt_einsum.contract as einsum
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsqr

import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(
  prog='Python Cheart Pfile Interface',
  description=
"""
Description to be added
""", formatter_class=RawTextHelpFormatter)
parser.add_argument('file', type=str,
  help='File name of data to be processed')

axes_dict = dict()
for i in range(8):
  for j in range(11):
    axes_dict[11*i + j] = (i,j)

def plot(file_name, time, data1, data2, label):
  fig, axs = plt.subplots(4, 1, dpi=300, figsize=(6, 8))
  axs[0].plot(time, data1[:,0,0], label="Old")
  axs[0].plot(time, data2[:,0,0], label="New")
  axs[0].set_title(fr'${label}_{{11}}$')
  axs[1].plot(time, data1[:,0,1])
  axs[1].plot(time, data2[:,0,1])
  axs[1].set_title(fr'${label}_{{12}}$')
  axs[2].plot(time, data1[:,1,0])
  axs[2].plot(time, data2[:,1,0])
  axs[2].set_title(fr'${label}_{{21}}$')
  axs[3].plot(time, data1[:,1,1])
  axs[3].plot(time, data2[:,1,1])
  axs[3].set_title(fr'${label}_{{22}}$')
  handles, labels = axs[0].get_legend_handles_labels()
  fig.legend(handles, labels, loc='upper right',prop={'size': 11})
  fig.tight_layout()
  plt.savefig(file_name, bbox_inches='tight')
  plt.close()


def plot_curve(fout, index, time, *data, label='S', unit='kPa', skip=1, msize=8):
  style = [
    {'marker': 'o', 'linestyle': 'None', 'mec': 'black', 'markersize': msize, 'markevery': skip,
     'fillstyle': 'none', 'markeredgewidth': .5, 'label' : " "},
    {'marker': 'None', 'linestyle': 'solid', 'color': 'red'},
    {'marker': 'None', 'linestyle': 'solid', 'color': 'blue'},
    {'marker': 'None', 'linestyle': 'solid', 'color': 'green'}]
  plt.rcParams['lines.linewidth'] = 1.0
  fig, axs = plt.subplots(8, 11, dpi=180, figsize=(44, 24))
  for k, (i,j) in enumerate(zip(index, index[1:])):
    for m, d in enumerate(data):
      axs[axes_dict[k]].plot(time[i:j], d[i:j], **style[m])
  fig.tight_layout()
  plt.savefig(fout, bbox_inches='tight')
  plt.close()


def get_n_data(index:ndarray, nelem=12, div=16):
  np = ascontiguousarray([(j-i) for (i,j) in zip(index, index[1:])])
  m  = stats.mode(np, keepdims=False).mode
  return ascontiguousarray([nelem if (2*n>m) else max(int(n/div),1) for n in np])


def make_joint_map(protocols:ndarray):
  n = 0
  map = list()
  for p in protocols:
    nodes = [max(0, n-1)]
    for _ in range(1, 2*(p+1)):
      n = n + 1
      nodes.append(n)
    map.append(ascontiguousarray(nodes))
  return map


def assemble_hermite_mat(time:ndarray, index, nelem=12, div=16):
  nelem = get_n_data(index, nelem, div)
  nnode = np.sum([(2*(n+1)) for n in nelem]) - (nelem.size - 1)
  map   = make_joint_map(nelem)
  A     = zeros((time.size, nnode), dtype=float)
  for k, (i, j) in enumerate(zip(index, index[1:])):
    A[i:j, map[k]] = get_hermite_matrix(x=time[i:j], t0=time[i], tend=time[min(j, index[-1]-1)], n=nelem[k])
  return nelem, map, csr_matrix(A)


class hermite_filter:

  def __init__(self, time:ndarray, index:ndarray, nelem=12, div=16) -> None:
    self.time  = time
    self.index = index
    self.end   = self.index[-1] - 1
    self.elems, self.map, self.A = assemble_hermite_mat(time, index, nelem, div)
    # self.Ainv  = linalg.pinv(A)
    self.intervals = list(zip(index, index[1:]))
    self.dof   = time.size - 2*self.elems.size
    self.all_vals = np.arange(0, index[-1])

  def interp(self, data:ndarray):
    # pars = einsum("ij,j->i", self.Ainv, data)
    pars = lsqr(self.A, data)[0]
    filtered = empty_like(data)
    for m, k, i, j in zip(self.map, self.elems, self.index, self.index[1:]):
      filtered[i:j] = hermite_interpolation(self.time[i:j], pars[m], self.time[i], self.time[min(j, self.end)], n=k)
    return filtered

  def filter(self, data:ndarray, eps_base:float=0.15):
    filtered = self.interp(data)
    eps = filtered - data
    eps = 2.1 * sqrt(einsum('i,i->', eps, eps) / self.dof) + eps_base
    bad_list = list()
    for i, j in zip(self.index, self.index[1:]):
      for m in range(i+1,j-1):
        if abs(data[m] - filtered[m]) > eps:
          for n in range(m-1, m+1):
            if (i < n < j-1):
              bad_list.append(n)
    nodes = np.setdiff1d(self.all_vals, bad_list)
    f = interpolate.interp1d(self.time[nodes], data[nodes])
    return f(self.time)

  def nested_filter(self, data:ndarray, eps_base:float=0.15, nest=5):
    fixed = data.copy()
    for _ in range(nest):
      fixed = self.filter(fixed, eps_base=eps_base)
    return fixed




def main(args=None):
  args = parser.parse_args(args=None)

  # ------------------------------------------------------------------------------------
  # Get data
  tstart = perf_counter()
  file_name = args.file
  tRead_1 = perf_counter()
  data  = pd.read_excel(file_name, sheet_name=None, index_col=None)
  tRead_2 = perf_counter()
  data  = pd.concat(data.values()).reset_index(drop=True)
  tRead_3 = perf_counter()
  print(f'Time to import excel:                   {tRead_2 - tRead_1:.7f}s')
  print(f'Time to convert join:                   {tRead_3 - tRead_2:.7f}s')
  nrows = data.shape[0]
  # ------------------------------------------------------------------------------------
  # Extract Data from Pandas
  tExtract_1 = perf_counter()
  time = data['Time_S'].to_numpy(dtype=float)
  x    = data[['X1', 'X2', 'X3', 'X4']].to_numpy(dtype=float)
  y    = data[['Y1', 'Y2', 'Y3', 'Y4']].to_numpy(dtype=float)
  coord= empty((nrows, 2, 4), dtype=float)
  coord[:,0,:] = x
  coord[:,1,:] = y
  f1   = data['XForce_mN'].to_numpy(dtype=float)
  f2   = data['YForce_mN'].to_numpy(dtype=float)
  # ------------------------------------------------------------------------------------
  # Extract indexes of different parts
  tExtract_2 = perf_counter()
  cycle_index = ascontiguousarray(np.append(data[['SetName', 'Cycle']].drop_duplicates().index.to_numpy(dtype=int), nrows))
  tExtract_3 = perf_counter()
  print(f'Time to data extraction:                {tExtract_2 - tExtract_1:.7f}s')
  print(f'Time to sorting index:                  {tExtract_3 - tExtract_2:.7f}s')


  tPInverse_0 = perf_counter()
  filter = hermite_filter(time=time, index=cycle_index, nelem=8)
  tPInverse_1 = perf_counter()

  print(f'{"Time to make the hermite filter:":<50}{tPInverse_1 - tPInverse_0:.7f}s')

  _, out_file = os.path.split(file_name)
  out_file, _ = os.path.splitext(out_file.replace(" ", "_"))
  os.makedirs(os.path.join('figures', out_file[:5]), exist_ok=True)

  for i in range(4):
    dat = x[:, i]
    fixed = filter.nested_filter(data=dat, eps_base=0.1, nest=6)
    filtered = filter.interp(fixed)
    plot_curve(os.path.join('figures', out_file[:5],
      f'filtered_x{i}.png'), cycle_index,
      time, dat, fixed, filtered, label='x')
    dat = y[:, i]
    fixed = filter.nested_filter(data=dat, eps_base=0.1, nest=6)
    filtered = filter.interp(fixed)
    plot_curve(os.path.join('figures', out_file[:5],
      f'filtered_y{i}.png'), cycle_index,
      time, dat, fixed, filtered, label='x')
  dat = f1[:]
  tFilter_0 = perf_counter()
  fixed = filter.nested_filter(data=dat, eps_base=2.5, nest=6)
  tFilter_1 = perf_counter()
  filtered = filter.interp(fixed)
  tFilter_2 = perf_counter()
  print(f'{"Time to filter data:":<50}{tFilter_1 - tFilter_0:.7f}s')
  print(f'{"Time to hermite interpolate data:":<50}{tFilter_2 - tFilter_1:.7f}s')
  plot_curve(os.path.join('figures', out_file[:5],
    f'filtered_f1.png'), cycle_index,
    time, dat, fixed, filtered, label='f')
  dat = f2[:]
  fixed = filter.nested_filter(data=dat, eps_base=2.5, nest=6)
  filtered = filter.interp(fixed)
  plot_curve(os.path.join('figures', out_file[:5],
    f'filtered_f2.png'), cycle_index,
    time, dat, fixed, filtered, label='f')

  tend = perf_counter()
  print(f'\nProgram completed in {tend - tstart} seconds.\n')


if __name__=="__main__":
  main()