# ruff: noqa: PLR0912, PLR0913, C901
import argparse
from argparse import RawTextHelpFormatter
from collections.abc import Callable, Sequence
from concurrent import futures
from math import sqrt
from pathlib import Path

import numpy as np
from numpy import matmul

from caputo_derivatives.analytical_solution import (
    diff_frac_p1_anal,
    diff_frac_p2_anal,
    diff_frac_t_anal,
    diff_t_anal,
    force_1,
    force_p,
    force_t,
    frac_t_anal,
)
from caputo_derivatives.caputo import CaputoInit, caputo_derivative1_array, caputo_diffeq1_array
from caputo_derivatives.tools.InputOutput import write_array
from caputo_derivatives.tools.ProgressBar import ProgressBar

parser = argparse.ArgumentParser(
    description="""
      To be made
      """,
    formatter_class=RawTextHelpFormatter,
)
parser.add_argument(
    "--clear-run",
    "-c",
    dest="cr",
    action="store_true",
    help="OPTIONAL: number of cores to compute with",
)
parser.add_argument(
    "--cores",
    "-n",
    dest="nc",
    type=int,
    default=1,
    help="OPTIONAL: number of cores to compute with",
)
_FArray1D = np.ndarray[tuple[int], np.dtype[np.float64]]
parm = np.array(
    [
        7.263437,
        -19.71696,
        100.7997,
        -509.0422,
        1106.24,
        -1204.629,
        718.3882,
        -236.194,
        39.39095,
        -2.497073,
    ]
)
parm2 = np.array(
    [
        6.38557,
        -8.97709,
        -12.97908,
        -77.38177,
        374.5206,
        -578.4844,
        449.7889,
        -191.9885,
        43.11304,
        -3.998423,
    ]
)


def test_func(
    func: Callable[[_FArray1D | Sequence[float], float, float, _FArray1D], _FArray1D],
    a: float,
    nk: int,
    d: float,
    tf: float,
    ti: int,
) -> _FArray1D:
    dt = 10.0 ** (-ti)
    nt = 2 * 10**ti + 1
    carp = CaputoInit(a, tf, nk)
    carp.delta = d
    time = np.arange(0.0, 2.0 + dt, dt)
    force = np.zeros(nt)
    if func.__name__ == "frac_t_anal":
        for j, t in enumerate(time):
            force[j] = force_t(parm, float(t))
        model = caputo_derivative1_array(force, dt, carp)
    elif func.__name__ == "diff_1_anal":
        for j, t in enumerate(time):
            force[j] = force_1(parm, float(t))
        model = caputo_diffeq1_array(force, dt, carp)
    elif func.__name__ == "diff_t_anal":
        for j, t in enumerate(time):
            force[j] = force_t(parm, float(t))
        model = caputo_diffeq1_array(force, dt, carp)
    elif func.__name__ == "diff_frac_t_anal":
        for j, t in enumerate(time):
            force[j] = force_t(parm, float(t))
        rhs = caputo_derivative1_array(force, dt, carp)
        model = caputo_diffeq1_array(rhs, dt, carp)
    elif func.__name__ == "diff_frac_p1_anal":
        for j, t in enumerate(time):
            force[j] = force_p(parm, float(t))
        rhs = caputo_derivative1_array(force, dt, carp)
        model = caputo_diffeq1_array(rhs, dt, carp)
    elif func.__name__ == "diff_frac_p2_anal":
        for j, t in enumerate(time):
            force[j] = force_p(parm2, float(t))
        rhs = caputo_derivative1_array(force, dt, carp)
        model = caputo_diffeq1_array(rhs, dt, carp)
    else:
        msg = "Function not recognized"
        raise ValueError(msg)
    return model


def run_test(
    func: Callable[[_FArray1D | Sequence[float], float, float, _FArray1D], _FArray1D],
    a: float,
    nk: int,
    d: float,
    tf: float,
    ti: int,
    sol: _FArray1D,
) -> float:
    dt = 10.0 ** (-ti)
    nt = 2 * 10**ti + 1
    model = np.zeros(nt)
    int_factor = np.ones(nt)
    int_factor[0] = 0.5
    int_factor[-1] = 0.5
    int_factor = dt * int_factor
    # Model Solution
    model = test_func(func, a, nk, d, tf, ti)
    res = model - sol
    res2 = res * res
    return sqrt(matmul(int_factor, res2))


def main(arg: list[str] | None = None) -> None:
    args = parser.parse_args(arg)
    with futures.ProcessPoolExecutor(args.nc) as executor:
        a_list = [0.1, 0.3, 0.5, 0.7, 0.9]
        n_list = [3, 6, 9, 12, 15]
        d_list = [0.05, 0.2, 0.8, 3.2, 12.8]
        e_list = [1, 2, 3, 4, 5]
        t_list = [2.0, 8.0, 32.0, 128.0, 512.0]
        tests = [frac_t_anal, diff_t_anal, diff_frac_t_anal, diff_frac_p1_anal, diff_frac_p2_anal]
        bar = ProgressBar(
            "Running Tf tests",
            len(a_list) * len(n_list) * len(d_list) * len(e_list) * len(t_list) * len(tests),
        )
        for fc in tests:
            for d in d_list:
                for a in a_list:
                    for nk in n_list:
                        future_to_time: dict[futures.Future[float], tuple[int, int]] = {}
                        fout = Path(
                            "data_tf", f"{fc.__name__}_{int(100 * d)}_{int(100 * a)}_{nk}.txt"
                        )
                        if not (args.cr) and fout.exists():
                            for _ in range(len(e_list) * len(t_list)):
                                bar.next()
                            continue
                        norm_arr = np.zeros((len(e_list), len(t_list)))
                        for i, ti in enumerate(e_list):
                            anal_sol = np.load(
                                Path(
                                    "analytical_solutions",
                                    f"{fc.__name__}_{int(100 * d)}_{int(100 * a)}_{ti}.dat",
                                )
                            )
                            for j, tf in enumerate(t_list):
                                future_to_time[
                                    executor.submit(run_test, fc, a, nk, d, tf, ti, anal_sol)
                                ] = (i, j)
                        for future in futures.as_completed(future_to_time):
                            i, j = future_to_time[future]
                            norm_arr[i, j] = future.result()
                            bar.next()
                        futures.wait(future_to_time, timeout=None)
                        write_array(fout, norm_arr)


if __name__ == "__main__":
    main()
