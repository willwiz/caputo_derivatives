import numpy as np
from numpy import exp

from .types import CaputoInit

type _ArrayCo[D: tuple[int, ...], T: (np.floating, np.integer)] = np.ndarray[D, np.dtype[T]]


def caputo_derivative1_array[F: np.floating, D: tuple[int, ...]](
    force: _ArrayCo[D, F], dt: float, carp: CaputoInit
) -> _ArrayCo[D, F]:
    dim = force.shape
    res = np.zeros_like(force)
    store = np.zeros(carp.n_p) if (len(dim) == 1) else np.zeros((carp.n_p, *dim[1:]))
    f_prev = force[0]
    c0 = carp.beta0 / dt
    coeff = np.zeros(carp.n_p)
    for k in range(carp.n_p):
        coeff[k] = carp.tau[k] / (carp.tau[k] + dt)
    for t in range(dim[0]):
        df = force[t] - f_prev
        f_prev = force[t]
        res[t] = c0 * df
        for k in range(carp.n_p):
            store[k] = coeff[k] * (store[k] + carp.beta[k] * df)
            res[t] = res[t] + store[k]
    return res


def caputo_derivative2_array[F: np.floating, D: tuple[int, ...]](
    force: _ArrayCo[D, F], dt: float, carp: CaputoInit
) -> _ArrayCo[D, F]:
    dim = force.shape
    res = np.zeros_like(force)
    store = np.zeros(carp.n_p) if (len(dim) == 1) else np.zeros((carp.n_p, *dim[1:]))
    f_prev = force[0]
    c0 = carp.beta0 / dt
    coeffq = np.zeros(carp.n_p)
    coeffb = np.zeros(carp.n_p)
    for k in range(carp.n_p):
        coeffb[k] = exp(-0.5 * dt / carp.tau[k])
        coeffq[k] = coeffb[k] * coeffb[k]
    for t in range(dim[0]):
        df = force[t] - f_prev
        f_prev = force[t]
        res[t] = c0 * df
        for k in range(carp.n_p):
            store[k] = coeffq[k] * store[k] + carp.beta[k] * coeffb[k] * df
            res[t] = res[t] + store[k]
    return res


def caputo_diffeq1_array[F: np.floating, D: tuple[int, ...]](
    force: _ArrayCo[D, F], dt: float, carp: CaputoInit
) -> _ArrayCo[D, F]:
    dim = force.shape
    res = np.zeros_like(force)
    store = np.zeros(carp.n_p) if (len(dim) == 1) else np.zeros((carp.n_p, *dim[1:]))
    f_prev = force[0]
    # Precompute
    c0 = carp.beta0 / dt
    coeff = np.zeros(carp.n_p)
    for k in range(carp.n_p):
        coeff[k] = carp.tau[k] / (carp.tau[k] + dt)
        c0 = c0 + carp.beta[k] * coeff[k]
    c0 = carp.delta * c0
    for t in range(dim[0]):
        res[t] = force[t] + c0 * f_prev
        for k in range(carp.n_p):
            res[t] = res[t] - carp.delta * coeff[k] * store[k]
        res[t] = res[t] / (1.0 + c0)
        # Update variables
        df = res[t] - f_prev
        f_prev = res[t]
        for k in range(carp.n_p):
            store[k] = coeff[k] * (store[k] + carp.beta[k] * df)
    return res


def caputo_diffeq2_array[F: np.floating, D: tuple[int, ...]](
    force: _ArrayCo[D, F], dt: float, carp: CaputoInit
) -> _ArrayCo[D, F]:
    dim = force.shape
    res = np.zeros_like(force)
    store = np.zeros(carp.n_p) if (len(dim) == 1) else np.zeros((carp.n_p, *dim[1:]))
    f_prev = force[0]
    # Precompute
    c0 = carp.beta0 / dt
    coeffb = np.zeros(carp.n_p)
    coeffq = np.zeros(carp.n_p)
    for k in range(carp.n_p):
        coeffb[k] = exp(-0.5 * dt / carp.tau[k])
        coeffq[k] = coeffb[k] * coeffb[k]
        c0 = c0 + carp.beta[k] * coeffb[k]
    c0 = carp.delta * c0
    for t in range(dim[0]):
        res[t] = force[t] + c0 * f_prev
        for k in range(carp.n_p):
            res[t] = res[t] - carp.delta * coeffq[k] * store[k]
        res[t] = res[t] / (1.0 + c0)
        # Update variables
        df = res[t] - f_prev
        f_prev = res[t]
        for k in range(carp.n_p):
            store[k] = coeffq[k] * store[k] + coeffb[k] * carp.beta[k] * df
    return res
