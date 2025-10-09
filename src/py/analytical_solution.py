from collections.abc import Sequence
from typing import TypeIs

import numpy as np
from scipy.special import gamma

from .mittag_leffler.mittag_leffler import ml

_FArray1D = np.ndarray[tuple[int], np.dtype[np.float64]]


def force_1(_a: _FArray1D, _t: float) -> float:
    return 1.0


def force_t(_a: _FArray1D, _t: float) -> float:
    return _t / 4


def force_p(a: _FArray1D, t: float) -> float:
    n = len(a)
    tk = 1.0
    res = 0.0
    for i in range(n):
        tk = tk * t
        res = res + a[i] * tk
    return res


def isfloat(x: object) -> TypeIs[float]:
    return isinstance(x, float)


def sametype[T: (_FArray1D, float)](_x: object, _y: T) -> TypeIs[T]:
    return True


# Bezare type narrowing from Pylance. WTF WHY?
def frac_t_anal[T: (_FArray1D, float)](
    _par: _FArray1D | Sequence[float], a: float, _d: float, t: T
) -> T:
    if isfloat(t):
        return t ** (1.0 - a) / gamma(2.0 - a) / 4.0
    return t ** (1.0 - a) / gamma(2.0 - a) / 4.0


def frac_tk_anal[T: (_FArray1D, float)](
    _par: _FArray1D | Sequence[float], a: float, k: float, t: T
) -> T:
    if isfloat(t):
        return t ** (k - a) * gamma(k + 1) / gamma(k + 1 - a)
    return t ** (k - a) * gamma(k + 1) / gamma(k + 1 - a)


def diff_1_anal[T: (_FArray1D, float)](
    _par: _FArray1D | Sequence[float], a: float, d: float, t: T
) -> T:
    if isfloat(t):
        g = t**a
        return g / d * ml(-g / d, a, a + 1.0)
    g = t**a
    return g / d * ml(-g / d, a, a + 1.0)


def diff_t_anal[T: (_FArray1D, float)](
    _par: _FArray1D | Sequence[float], a: float, d: float, t: T
) -> T:
    if isfloat(t):
        return t ** (a + 1.0) / d * ml(-(t**a) / d, a, a + 2.0) / 4
    return t ** (a + 1.0) / d * ml(-(t**a) / d, a, a + 2.0) / 4


def diff_frac_t_anal[T: (_FArray1D, float)](
    _par: _FArray1D | Sequence[float], a: float, d: float, t: T
) -> T:
    if isfloat(t):
        return t / d * ml(-(t**a) / d, a, 2.0) / 4
    return t / d * ml(-(t**a) / d, a, 2.0) / 4


def diff_frac_p_anal[T: (_FArray1D, float)](
    par: _FArray1D | Sequence[float], a: float, d: float, t: T
) -> T:
    tk = 1.0
    res = 0.0 * t
    for i in range(len(par)):
        tk = tk * t
        res = res + par[i] * tk * gamma(2.0 + i) / d * ml(-(t**a) / d, a, i + 2.0)
    if sametype(res, t):
        return res
    msg = "It's insane that I have to do this"
    raise TypeError(msg)


def diff_frac_p1_anal[T: (_FArray1D, float)](
    par: _FArray1D | Sequence[float], a: float, d: float, t: T
) -> T:
    tk = 1.0
    res = 0.0
    for i in range(len(par)):
        tk = tk * t
        res = res + par[i] * tk * gamma(2.0 + i) / d * ml(-(t**a) / d, a, i + 2.0)
    if sametype(res, t):
        return res
    msg = "It's insane that I have to do this"
    raise TypeError(msg)


def diff_frac_p2_anal[T: (_FArray1D, float)](
    par: _FArray1D | Sequence[float], a: float, d: float, t: T
) -> T:
    tk = 1.0
    res = 0.0
    for i in range(len(par)):
        tk = tk * t
        res = res + par[i] * tk * gamma(2.0 + i) / d * ml(-(t**a) / d, a, i + 2.0)
    if sametype(res, t):
        return res
    msg = "It's insane that I have to do this"
    raise TypeError(msg)


def diff_frac_p3_anal[T: (_FArray1D, float)](
    par: _FArray1D | Sequence[float], a: float, d: float, t: T
) -> T:
    tk = 1.0
    res = 0.0
    for i in range(len(par)):
        tk = tk * t
        res = res + par[i] * tk * gamma(2.0 + i) / d * ml(-(t**a) / d, a, i + 2.0)
    if sametype(res, t):
        return res
    msg = "It's insane that I have to do this"
    raise TypeError(msg)
