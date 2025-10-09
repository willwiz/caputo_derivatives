import dataclasses
import unittest

import numpy as np
from numpy.testing import assert_almost_equal

from py.matlaws.holzapfel_models import (
    NeoHookean2D,
)
from src.py.matlaws.neohookean import NeoHookeanPy

# from src.py.caputo import

np.set_printoptions(suppress=True)
__unittest = True


@dataclasses.dataclass
class Benchmark:
    time: np.ndarray
    strain: np.ndarray
    neohookean: np.ndarray
    hog: np.ndarray
    hogstruc: np.ndarray


class ResidualCalculations(unittest.TestCase):
    def __init__(self) -> None:
        super().__init__()
        time = np.linspace(0, 1.0, 11)
        identity = np.array([1.0, 0.0, 0.0, 1.0])
        max_strain = np.array([0.2, 0.0, 0.0, 0.4])
        self.strain = np.einsum("i,j->ij", time, max_strain) + identity

    def test_hyperelastic(self) -> None:
        neo_cpp = np.empty_like(self.strain)
        neo_py = np.empty_like(self.strain)
        model_py = NeoHookeanPy(0.5)
        for i, v in enumerate(self.strain):
            neo_py[i] = model_py.stress(v)
        model_cpp = NeoHookean2D(0.5)
        for i, v in enumerate(self.strain):
            neo_cpp[i] = model_cpp.stress(v)
        assert_almost_equal(neo_cpp, neo_py, decimal=10)


if __name__ == "__main__":
    pass
