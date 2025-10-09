import numpy as np

from py.caputo._io import import_precomputed_data

_FArray1D = np.ndarray[tuple[int], np.dtype[np.float64]]


class CaputoInit:
    __slots__ = ("alpha", "beta", "beta0", "delta", "freq", "n_p", "tau")
    beta: _FArray1D
    tau: _FArray1D
    beta0: float
    freq: float
    n_p: int

    def __init__(self, alpha: float, tf: float, n_p: int) -> None:
        freq = 2.0 * np.pi / tf
        carp = import_precomputed_data(n_p)
        betam = carp["betam"]
        tau = carp["taum"]
        beta: np.ndarray[tuple[int, int], np.dtype[np.float64]] = betam[:-1]
        x = np.linspace(0.01, 1.0, num=100, endpoint=True)
        self.n_p = n_p
        self.freq = freq
        self.beta0 = np.interp(alpha, x, betam[-1]) * (freq ** (alpha - 1.0))
        self.beta = np.array([np.interp(alpha, x, i) * (freq**alpha) for i in beta])
        self.tau = np.array([np.interp(alpha, x, i) / freq for i in tau])
