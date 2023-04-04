import unittest
__unittest = True

import numpy as np
import numpy as np
import numpy.typing as npt
from numpy.testing import assert_almost_equal
import typing as tp
from time import perf_counter


from src.py.caputoD import *
from src.py.matlaws.HModels import (
    caputo_initialize,
)
from src.py.AnalyticalSolution import (
    frac_tk_anal,
    frac_t_anal,
    diff_frac_p_anal,
)


np.set_printoptions(suppress=True)


def frac_deriv_cpp(alpha:float, dt:float, data:npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    res = np.zeros_like(data)
    n = len(data)
    carp = caputo_initialize(alpha, 10.0, 0.0)
    for i in range(1, n):
        res[i] = carp.caputo_iter(data[i], dt)
    return res

def diffeq_deriv_cpp(alpha:float, delta:float, dt:float, data:npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    rhs = np.zeros_like(data)
    res = np.zeros_like(data)
    n = len(data)
    carpR = caputo_initialize(alpha, 10.0, delta)
    carpL = caputo_initialize(alpha, 10.0, delta)
    for i in range(1, n):
        rhs[i] = carpR.caputo_iter(data[i], dt)
        res[i] = carpL.diffeq_iter(rhs[i], dt)
    return res

def frac_deriv_py(alpha:float, dt:float, data:npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    res = np.zeros_like(data)
    carp = caput_init(alpha, 10.0, 9)
    res = caputo_derivative2_array(data, dt, carp)
    return res

def diffeq_deriv_py(alpha:float, delta:float, dt:float, data:npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    res = np.zeros_like(data)
    carpR = caput_init(alpha, 10.0, 9)
    carpL = caput_init(alpha, 10.0, 9)
    carpL.delta = delta
    rhs = caputo_derivative2_array(data, dt, carpR)
    res = caputo_diffeq2_array(rhs, dt, carpL)
    return res

# ------------------------------------------------------------------------------
# Unit test class
# ------------------------------------------------------------------------------
class hyperelastic_simulations(unittest.TestCase):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        file_name = 'unittests/benchmark_caputo'
        datazip = np.load(f'{file_name}.npz')
        self.alpha  = datazip['alpha']
        self.delta  = datazip['delta']
        self.dt     = datazip['dt']
        self.time   = datazip['time']
        self.x      = datazip['x']
        self.y      = datazip['y']
        self.diffeq = datazip['diffeq']
    # The Unit tests
    def test_cvp(self):
        pyt = frac_deriv_py(self.alpha, self.dt, self.x)
        cpp = frac_deriv_cpp(self.alpha, self.dt, self.x)
        assert_almost_equal(cpp, pyt)
    def test_py(self):
        cpp = frac_deriv_py(self.alpha, self.dt, self.x)
        norm = np.linalg.norm(cpp-self.y, ord='fro')
        self.assertAlmostEqual(norm, 0.0, 2)
    def test_cpp(self):
        cpp = frac_deriv_cpp(self.alpha, self.dt, self.x)
        norm = np.linalg.norm(cpp-self.y, ord='fro')
        self.assertAlmostEqual(norm, 0.0, 2)

    def test_cvp_diffeq(self):
        pyt = diffeq_deriv_py(self.alpha, self.delta, self.dt, self.x)
        cpp = diffeq_deriv_cpp(self.alpha, self.delta, self.dt, self.x)
        assert_almost_equal(cpp, pyt)
    def test_py_diffeq(self):
        cpp = diffeq_deriv_py(self.alpha, self.delta, self.dt, self.x)
        norm = np.linalg.norm(cpp - self.diffeq, ord='fro')
        self.assertAlmostEqual(norm, 0.0, 2)
    def test_cpp_diffeq(self):
        cpp = diffeq_deriv_cpp(self.alpha, self.delta, self.dt, self.x)
        norm = np.linalg.norm(cpp - self.diffeq, ord='fro')
        self.assertAlmostEqual(norm, 0.0, 2)

def gen_testdata():
    dt = 0.001
    nt = round(1/dt)+1
    alpha = 0.2
    delta = 0.05
    time = np.linspace(0, 1.0, round(1/dt)+1)
    data = np.zeros((nt,4), dtype=float)
    res = np.zeros_like(data)
    diffeq = np.zeros_like(data)
    for i in range(4):
        data[:, i] = time**(i+1)
        res[:, i] = frac_tk_anal(0.0, alpha, i+1, time)
    diffeq[:, 0] = diff_frac_p_anal([1,0,0,0], alpha, delta, time)
    diffeq[:, 1] = diff_frac_p_anal([0,1,0,0], alpha, delta, time)
    diffeq[:, 2] = diff_frac_p_anal([0,0,1,0], alpha, delta, time)
    diffeq[:, 3] = diff_frac_p_anal([0,0,0,1], alpha, delta, time)
    np.savez('unittests/benchmark_caputo', alpha=alpha, delta = delta, dt=dt, time=time, x = data, y = res, diffeq=diffeq)
    return



if __name__=="__main__":
    gen_testdata()