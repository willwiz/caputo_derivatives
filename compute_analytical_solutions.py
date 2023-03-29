#!/usr/bin/env python3

import numpy as np
from src.py.caputoD import *
from src.py.AnalyticalSolution import *
from src.py.tools.ProgressBar import progress_bar
from os.path import exists
from typing import Callable
from concurrent import futures

import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description="""
      To be made
      """, formatter_class=RawTextHelpFormatter)
parser.add_argument('--n-cores',              '-n',   dest='nc',  type=int, default=1,
                    help='OPTIONAL: number of cores to compute with')


parm  = np.array([7.263437, -19.71696, 100.7997, -509.0422, 1106.24, -1204.629, 718.3882, -236.194, 39.39095, -2.497073])
parm2 = np.array([6.38557, -8.97709, -12.97908, -77.38177, 374.5206, -578.4844, 449.7889, -191.9885, 43.11304, -3.998423])


def calc_sol(exec, func: Callable, par : np.ndarray, a : float, d : float) -> None:
  e_list = [1,2,3,4,5]
  for i, ti in enumerate(e_list):
    dt = 10.0**(-ti)
    nt = 2*10**ti + 1
    fout = "analytical_solutions/"+f"{func.__name__}_{int(100*a)}_{int(100*d)}_{i+1}.dat"
    if exists(fout) : continue
    time = np.arange(0.0, 2.0+dt, dt)
    anal_sol = np.zeros(nt)
    bar = progress_bar(f"precomputing analytical tests for a = {a}, d = {d}, dt = {dt}", len(time))
    future_to_time = {exec.submit(func, par, a, d, time[k]): k for k in range(len(time))}
    for future in futures.as_completed(future_to_time):
      k = future_to_time[future]
      anal_sol[k] = future.result()
      bar.next()
    with open(fout, 'wb') as f:
      np.save(f, anal_sol)

def submit_calc(exec, func: Callable, par : np.ndarray, a : float, d : float, i : int, fout) -> None:
  dt = 10.0**(-i)
  nt = 2*10**i + 1
  time = np.arange(0.0, 2.0+dt, dt)
  anal_sol = np.zeros(nt)
  bar = progress_bar(f"{func.__name__:<17} with d = {d:<4}, a = {a:<3}, dt = {dt:<6}", len(time))
  future_to_time = {exec.submit(func, par, a, d, time[k]): k for k in range(len(time))}
  for future in futures.as_completed(future_to_time):
    k = future_to_time[future]
    anal_sol[k] = future.result()
    bar.next()
  with open(fout, 'wb') as f:
    np.save(f, anal_sol)


def main(args=None):
  args=parser.parse_args(args)
  with futures.ProcessPoolExecutor(args.nc) as executor:
    a_list  = [0.1, 0.3, 0.5, 0.7, 0.9]
    d_list  = [0.05, 0.2, 0.8, 3.2, 12.8]
    e_list  = [1, 2, 3, 4, 5]
    tests   = [
      [frac_t_anal, parm],
      [diff_t_anal, parm],
      [diff_frac_t_anal, parm],
      [diff_frac_p1_anal, parm],
      [diff_frac_p2_anal, parm2]
    ]
    for fc, par in tests:
      for d in d_list:
        for a in a_list:
          for i in e_list:
            fout = "analytical_solutions/"+f"{fc.__name__}_{int(100*d)}_{int(100*a)}_{i}.dat"
            if exists(fout) : continue
            submit_calc(executor, fc, par, a, d, i, fout)

if __name__=="__main__":
  main()