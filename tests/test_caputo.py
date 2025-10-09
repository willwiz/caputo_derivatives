# ruff: noqa: S101
import unittest

import numpy as np
from numpy.testing import assert_almost_equal

from py.analytical_solution import (
    diff_frac_p_anal,
    frac_tk_anal,
)
from py.matlaws.holzapfel_models import (
    caputo_initialize,
)
from src.py.caputo import CaputoInit, caputo_derivative2_array, caputo_diffeq2_array

np.set_printoptions(suppress=True)

__unittest = True
_FArray1D = np.ndarray[tuple[int], np.dtype[np.float64]]


def frac_deriv_cpp(alpha: float, dt: float, data: _FArray1D) -> _FArray1D:
    res = np.zeros_like(data)
    n = len(data)
    carp = caputo_initialize(alpha, 10.0, 0.0)
    for i in range(1, n):
        res[i] = carp.caputo_iter(data[i], dt)
    return res


def diffeq_deriv_cpp(alpha: float, delta: float, dt: float, data: _FArray1D) -> _FArray1D:
    rhs = np.zeros_like(data)
    res = np.zeros_like(data)
    n = len(data)
    carp_r = caputo_initialize(alpha, 10.0, delta)
    carp_l = caputo_initialize(alpha, 10.0, delta)
    for i in range(1, n):
        rhs[i] = carp_r.caputo_iter(data[i], dt)
        res[i] = carp_l.diffeq_iter(rhs[i], dt)
    return res


def frac_deriv_py(alpha: float, dt: float, data: _FArray1D) -> _FArray1D:
    carp = CaputoInit(alpha, 10.0, 9)
    return caputo_derivative2_array(data, dt, carp)


def diffeq_deriv_py(alpha: float, delta: float, dt: float, data: _FArray1D) -> _FArray1D:
    carp_r = CaputoInit(alpha, 10.0, 9)
    carp_l = CaputoInit(alpha, 10.0, 9)
    carp_l.delta = delta
    rhs = caputo_derivative2_array(data, dt, carp_r)
    return caputo_diffeq2_array(rhs, dt, carp_l)


# ------------------------------------------------------------------------------
# Unit test class
# ------------------------------------------------------------------------------
TOL = 1e-2


class HyperelasticSimulations(unittest.TestCase):
    def __init__(self) -> None:
        super().__init__()
        file_name = "unittests/benchmark_caputo"
        datazip = np.load(f"{file_name}.npz")
        self.alpha = datazip["alpha"]
        self.delta = datazip["delta"]
        self.dt = datazip["dt"]
        self.time = datazip["time"]
        self.x = datazip["x"]
        self.y = datazip["y"]
        self.diffeq = datazip["diffeq"]

    # The Unit tests
    def test_cvp(self) -> None:
        pyt = frac_deriv_py(self.alpha, self.dt, self.x)
        cpp = frac_deriv_cpp(self.alpha, self.dt, self.x)
        assert_almost_equal(cpp, pyt)

    def test_py(self) -> None:
        cpp = frac_deriv_py(self.alpha, self.dt, self.x)
        norm = np.linalg.norm(cpp - self.y, ord="fro")
        assert abs(norm) < TOL

    def test_cpp(self) -> None:
        cpp = frac_deriv_cpp(self.alpha, self.dt, self.x)
        norm = np.linalg.norm(cpp - self.y, ord="fro")
        assert abs(norm) < TOL

    def test_cvp_diffeq(self) -> None:
        pyt = diffeq_deriv_py(self.alpha, self.delta, self.dt, self.x)
        cpp = diffeq_deriv_cpp(self.alpha, self.delta, self.dt, self.x)
        assert_almost_equal(cpp, pyt)

    def test_py_diffeq(self) -> None:
        cpp = diffeq_deriv_py(self.alpha, self.delta, self.dt, self.x)
        norm = np.linalg.norm(cpp - self.diffeq, ord="fro")
        assert abs(norm) < TOL

    def test_cpp_diffeq(self) -> None:
        cpp = diffeq_deriv_cpp(self.alpha, self.delta, self.dt, self.x)
        norm = np.linalg.norm(cpp - self.diffeq, ord="fro")
        assert abs(norm) < TOL


def gen_testdata() -> None:
    dt = 0.001
    nt = round(1 / dt) + 1
    alpha = 0.2
    delta = 0.05
    time = np.linspace(0, 1.0, round(1 / dt) + 1)
    data = np.zeros((nt, 4), dtype=float)
    res = np.zeros_like(data)
    diffeq = np.zeros_like(data)
    for i in range(4):
        data[:, i] = time ** (i + 1)
        res[:, i] = frac_tk_anal([0.0], alpha, i + 1, time)
    diffeq[:, 0] = diff_frac_p_anal([1, 0, 0, 0], alpha, delta, time)
    diffeq[:, 1] = diff_frac_p_anal([0, 1, 0, 0], alpha, delta, time)
    diffeq[:, 2] = diff_frac_p_anal([0, 0, 1, 0], alpha, delta, time)
    diffeq[:, 3] = diff_frac_p_anal([0, 0, 0, 1], alpha, delta, time)
    np.savez(
        "unittests/benchmark_caputo",
        alpha=alpha,
        delta=delta,
        dt=dt,
        time=time,
        x=data,
        y=res,
        diffeq=diffeq,
    )


if __name__ == "__main__":
    gen_testdata()
