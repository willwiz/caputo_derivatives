import numpy as np
from scipy.special import gamma

from .ml_internal import LTInversion


def ml(z, alpha, beta=1.0, gama=1.0):
    eps = np.finfo(np.float64).eps
    if (
        np.real(alpha) <= 0
        or np.real(gama) <= 0
        or np.imag(alpha) != 0.0
        or np.imag(beta) != 0.0
        or np.imag(gama) != 0.0
    ):
        msg = "ALPHA and GAMA must be real and positive. BETA must be real."
        raise ValueError(msg)
    if np.abs(gama - 1.0) > eps:
        if alpha > 1.0:
            msg = "GAMMA != 1 requires 0 < ALPHA < 1"
            raise ValueError(msg)
        if (np.abs(np.angle(np.repeat(z, np.abs(z) > eps))) <= alpha * np.pi).any():
            msg = "|Arg(z)| <= alpha*pi"
            raise ValueError(msg)

    return np.vectorize(ml_, [np.float64])(z, alpha, beta, gama)


def ml_(z, alpha, beta, gama):
    # Target precision
    log_epsilon = np.log(1.0e-15)
    # Inversion of the LT
    if np.abs(z) < 1.0e-15:
        return 1 / gamma(beta)
    return LTInversion(1, z, alpha, beta, gama, log_epsilon)
