# ruff: noqa: N801, N802
import numpy as np

__all__ = [
    "CaputoSimulation_C_M",
    "CaputoSimulation_C_M_Scaled",
    "GetScaledParametersScaled_cpp",
    "GetScaledParameters_cpp",
    "HESimulation",
    "HESimulation_Scaled",
    "HOGstruc2D",
    "NeoHookean2D",
    "caputo_initialize",
    "pyHOG2D",
]
type _FArray1D = np.ndarray[tuple[int], np.dtype[np.float64]]
type _FArray2D = np.ndarray[tuple[int, int], np.dtype[np.float64]]

class caputo_initialize:
    def __init__(self, alpha: float, tf: float, delta: float) -> None: ...
    def set_pars(self, alpha: float, tf: float, delta: float) -> None: ...
    @property
    def beta0(self) -> float: ...
    @property
    def betas(self) -> _FArray1D: ...
    @property
    def taus(self) -> _FArray1D: ...
    def caputo_iter(self, f: _FArray1D, dt: float) -> _FArray1D: ...
    def diffeq_iter(self, f: _FArray1D, dt: float) -> _FArray1D: ...

class NeoHookean2D:
    def __init__(self, mu: float) -> None: ...
    def set_pars(self, mu: float) -> None: ...
    def stress(self, f: _FArray1D) -> _FArray1D: ...

class HOGstruc2D:
    def __init__(
        self, k1: float, k2: float, theta: float, alpha: float, beta: float, kip: float, kop: float
    ) -> None: ...
    def set_pars(
        self, k1: float, k2: float, theta: float, alpha: float, beta: float, kip: float, kop: float
    ) -> None: ...
    def stress(self, f: _FArray1D) -> _FArray1D: ...

class pyHOG2D:
    def __init__(self, k1: float, k2: float, theta: float) -> None: ...
    def set_pars(self, k1: float, k2: float, theta: float) -> None: ...
    def stress(self, f: _FArray1D) -> _FArray1D: ...

def GetScaledParameters_cpp(
    pars: _FArray1D, fiber: _FArray1D, visco: _FArray1D, tf: float
) -> _FArray1D: ...
def GetScaledParametersScaled_cpp(
    pars: _FArray1D, fiber: _FArray1D, visco: _FArray1D, tf: float, cmax: _FArray1D
) -> _FArray1D: ...
def HESimulation(
    pars: _FArray1D,
    fiber: _FArray1D,
    caputo: _FArray1D,
    tf: float,
    strain: _FArray2D,
    dt: _FArray1D,
) -> _FArray2D: ...
def HESimulation_Scaled(
    pars: _FArray1D,
    fiber: _FArray1D,
    caputo: _FArray1D,
    tf: float,
    cmax: _FArray1D,
    strain: _FArray2D,
    dt: _FArray1D,
) -> _FArray2D: ...
def CaputoSimulation_C_M(
    pars: _FArray1D,
    fiber: _FArray1D,
    visco: _FArray1D,
    tf: float,
    strain: _FArray2D,
    dt: _FArray1D,
) -> _FArray2D: ...
def CaputoSimulation_C_M_Scaled(
    pars: _FArray1D,
    fiber: _FArray1D,
    visco: _FArray1D,
    tf: float,
    cmax: _FArray1D,
    strain: _FArray2D,
    dt: _FArray1D,
) -> _FArray2D: ...
