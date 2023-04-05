import dataclasses
import numpy as np
from numpy import zeros, array
from numpy.testing import assert_almost_equal
from src.py.caputoD import *
from src.py.matlaws.neohookean import (
    NeoHookean_py
)
from src.py.matlaws.HModels import (
    NeoHookean2D,
)
from time import perf_counter
import unittest

np.set_printoptions(suppress=True)
__unittest = True

@dataclasses.dataclass
class benchmark:
    time:np.ndarray
    strain:np.ndarray
    neohookean:np.ndarray
    hog:np.ndarray
    hogstruc:np.ndarray


class residual_calculations(unittest.TestCase):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        time = np.linspace(0, 1.0, 11)
        identity = np.array([1.0,0.0,0.0,1.0])
        max_strain = np.array([0.2,0.0,0.0,0.4])
        self.strain = np.einsum('i,j->ij', time, max_strain) + identity
    def test_Hyperelastic(self):
        neo_cpp = np.empty_like(self.strain)
        neo_py = np.empty_like(self.strain)
        model_py = NeoHookean_py(0.5)
        for i, v in enumerate(self.strain):
            neo_py[i] = model_py.stress(v)
        model_cpp = NeoHookean2D(0.5)
        for i, v in enumerate(self.strain):
            neo_cpp[i] = model_cpp.stress(v)
        assert_almost_equal(neo_cpp, neo_py, decimal=10)


if __name__=="__main__":
    pass