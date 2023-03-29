#!/usr/bin/env python3

import numpy as np
from numpy import matmul
from math import sqrt
from src.py.caputoD import *
from src.py.AnalyticalSolution import *
from src.py.tools.ProgressBar import progress_bar
from src.py.tools.InputOutput import write_array
from concurrent import futures
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description="""
      To be made
      """, formatter_class=RawTextHelpFormatter)
parser.add_argument('--clear-run',              '-c',   dest='cr', action='store_true',
                    help='OPTIONAL: number of cores to compute with')
parser.add_argument('--cores',              '-n',   dest='nc',  type=int, default=1,
                    help='OPTIONAL: number of cores to compute with')

parm  = np.array([7.263437, -19.71696, 100.7997, -509.0422, 1106.24, -1204.629, 718.3882, -236.194, 39.39095, -2.497073])
parm2 = np.array([6.38557, -8.97709, -12.97908, -77.38177, 374.5206, -578.4844, 449.7889, -191.9885, 43.11304, -3.998423])


def test_func(func, a, nk, d, tf, ti):
  dt = 10.0**(-ti)
  nt = 2*10**ti + 1
  carp = caput_init(a, tf, nk)
  carp.delta = d
  time = np.arange(0.0, 2.0+dt, dt)
  force    = np.zeros(nt)
  if func.__name__ == 'frac_t_anal':
    for j, t in enumerate(time):
      force[j] = force_t(parm, t)
    model = caputo_derivative1_array(force, dt, carp)
  elif func.__name__ == 'diff_1_anal':
    for j, t in enumerate(time):
      force[j] = force_1(parm, t)
    model = caputo_diffeq1_array(force, dt, carp)
  elif func.__name__ == 'diff_t_anal':
    for j, t in enumerate(time):
      force[j] = force_t(parm, t)
    model = caputo_diffeq1_array(force, dt, carp)
  elif func.__name__ == 'diff_frac_t_anal':
    for j, t in enumerate(time):
      force[j] = force_t(parm, t)
    rhs   = caputo_derivative1_array(force, dt, carp)
    model = caputo_diffeq1_array(rhs, dt, carp)
  elif func.__name__ == 'diff_frac_p1_anal':
    for j, t in enumerate(time):
      force[j] = force_p(parm, t)
    rhs   = caputo_derivative1_array(force, dt, carp)
    model = caputo_diffeq1_array(rhs, dt, carp)
  elif func.__name__ == 'diff_frac_p2_anal':
    for j, t in enumerate(time):
      force[j] = force_p(parm2, t)
    rhs   = caputo_derivative1_array(force, dt, carp)
    model = caputo_diffeq1_array(rhs, dt, carp)
  return model


def run_test(func, a, nk, d, tf, ti, sol):
  dt = 10.0**(-ti)
  nt = 2*10**ti + 1
  model    = np.zeros(nt)
  int_factor = np.ones(nt)
  int_factor[0]  = 0.5
  int_factor[-1] = 0.5
  int_factor = dt * int_factor
  # Model Solution
  model = test_func(func, a, nk, d, tf, ti)
  res = model - sol
  res2 = res * res
  return sqrt( matmul(int_factor, res2))


def main(args=None):
  args=parser.parse_args(args)
  with futures.ProcessPoolExecutor(args.nc) as executor:
    a_list  = [0.1, 0.3, 0.5, 0.7, 0.9]
    n_list  = [3, 6, 9, 12, 15]
    d_list  = [0.05, 0.2, 0.8, 3.2, 12.8]
    e_list  = [1, 2, 3, 4, 5]
    t_list  = [32.0]
    tests   = [
      frac_t_anal,
      diff_t_anal,
      diff_frac_t_anal,
      diff_frac_p1_anal,
      diff_frac_p2_anal
      ]
    bar = progress_bar("Running Prony tests", len(a_list)*len(n_list)*len(d_list)*len(e_list)*len(t_list)*len(tests))
    for fc in tests:
      for d in d_list:
        for a in a_list:
          for tf in t_list:
            future_to_time = {}
            fout = os.path.join("data_prony", f"{fc.__name__}_{int(100*d)}_{int(100*a)}_{int(tf)}.txt")
            if not(args.cr) and os.path.exists(fout):
              for _ in range(len(n_list)*len(e_list)):
                bar.next()
              continue
            norm_arr = np.zeros((len(n_list),len(e_list)))
            for i, ti in enumerate(e_list):
              anal_sol = np.load(os.path.join("analytical_solutions", f'{fc.__name__}_{int(100*d)}_{int(100*a)}_{ti}.dat'))
              for j, nk in enumerate(n_list):
                future_to_time[executor.submit(run_test, fc, a, nk, d, tf, ti, anal_sol)] = (i,j)
            for future in futures.as_completed(future_to_time):
              i, j = future_to_time[future]
              norm_arr[j, i] = future.result()
              bar.next()
            futures.wait(future_to_time, timeout=None)
            write_array(fout, norm_arr)
  return

if __name__=="__main__":
  main()