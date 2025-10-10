# File: model.pyx
# ruff: noqa: N806
"""Implements the Holzapfel Gasser Ogden model 2013 version.

Updated in 2015 in Holzapfel et al. Modelling non-symmetric collagen ...
J. Royal Soc. Interface. 2015.

Jan 15, 2019
By: Will Zhang

HOG-2013 Model
Parameter List
mu      :   Matrix modulus
k1      :   Collagen modulus
k2      :   Collagen exponent
e1      :   Elastin modulus
e2      :   Elastin exponent
s1      :   Smooth muscle exponent
s2      :   Smooth muscle exponent

Vector
alpha   :   Fiber split
kop     :   splay parameters
kip     :   splay parameters

Inputs:
C11
C12
C21
C22

Outputs:
S11
S12
S21
S22

"""

from collections.abc import Sequence

import numpy as np
from numpy import array, cos, exp, ndarray, sin

_FArray1D = ndarray[tuple[int], np.dtype[np.float64]]

Identity2D: _FArray1D = array([[1, 0, 0, 1]], dtype=float)


class Hog20132D:
    __slots__ = (
        "He",
        "Hs",
        "a",
        "b",
        "c",
        "e1",
        "e2",
        "h4",
        "h6",
        "k1",
        "k2",
        "m4",
        "m6",
        "mu",
        "s1",
        "s2",
    )

    mu: float
    k1: float
    k2: float
    e1: float
    e2: float
    s1: float
    s2: float
    a: float
    b: float
    c: float
    h4: _FArray1D
    h6: _FArray1D
    He: _FArray1D
    Hs: _FArray1D

    def __init__(self, pars: _FArray1D, alpha: float, kip: float, kop: float) -> None:
        # unpack parameters
        self.mu = pars[0]
        self.k1 = pars[1]
        self.k2 = pars[2]
        self.e1 = pars[3]
        self.e2 = pars[4]
        self.s1 = pars[5]
        self.s2 = pars[6]
        # Compute fiber info
        self.a = 2.0 * kop * kip
        self.b = 2.0 * kop * (1.0 - 2.0 * kip)
        self.c = 1.0 - 3.0 * self.a - self.b
        # Compute the othogonal fiber vector
        cosa = cos(alpha)
        sina = sin(alpha)
        self.m4 = array([cosa * cosa, cosa * sina, cosa * sina, sina * sina])
        self.m6 = array([cosa * cosa, -cosa * sina, -cosa * sina, sina * sina])
        self.h4 = self.a * Identity2D + self.b * self.m4
        self.h6 = self.a * Identity2D + self.b * self.m6
        # print(self.H4 + self.H6)
        self.He = array([1, 0, 0, 0], dtype=float)
        self.Hs = array([0, 0, 0, 1], dtype=float)

    def stress(self, args: _FArray1D) -> _FArray1D:
        det: float = args[0] * args[3] - args[1] * args[1]
        i_n = 1.0 / det
        c_inv: _FArray1D = array([args[3], -args[1], -args[2], args[0]]) * i_n
        I_1 = args[0] + args[3] + i_n
        I_4 = self.m4 @ args
        I_6 = self.m6 @ args
        e6 = self.a * I_1 + self.c * i_n - 1.0
        e4 = e6 + self.b * I_4
        e6 = e6 + self.b * I_6
        dWd4 = self.k1 * e4 * exp(self.k2 * e4 * e4)
        dWd6 = self.k1 * e6 * exp(self.k2 * e6 * e6)
        # print(I_4, I_6, dWd4, dWd6)
        Ee = args[0] - 1.0
        Es = args[3] - 1.0
        dWde = self.e1 * Ee * exp(self.e2 * Ee * Ee)
        dWds = self.s1 * Es * exp(self.s2 * Es * Es)
        p = (self.mu + self.c * dWd4 + self.c * dWd6) * i_n
        val = (
            self.mu * Identity2D + dWd4 * self.h4 + dWd6 * self.h6 + dWde * self.He + dWds * self.Hs
        )
        return val - p * c_inv

    def stress_parts(self, args: _FArray1D) -> Sequence[_FArray1D]:
        det = args[0] * args[3] - args[1] * args[1]
        I_n = 1.0 / det
        Cinv = array([args[3], -args[1], -args[2], args[0]]) * I_n
        I_1 = args[0] + args[3] + I_n
        I_4 = self.m4 @ args
        I_6 = self.m6 @ args
        e6 = self.a * I_1 + self.c * I_n - 1.0
        e4 = e6 + self.b * I_4
        e6 = e6 + self.b * I_6
        dWd4 = self.k1 * e4 * exp(self.k2 * e4 * e4)
        dWd6 = self.k1 * e6 * exp(self.k2 * e6 * e6)
        # print(I_4, I_6, dWd4, dWd6)
        Ee = args[0] - 1.0
        Es = args[3] - 1.0
        dWde = self.e1 * Ee * exp(self.e2 * Ee * Ee)
        dWds = self.s1 * Es * exp(self.s2 * Es * Es)
        # p = (self.mu + self.C * dWd4 + self.C * dWd6)
        # val = self.mu * identity2D + dWd4*self.H4 + dWd6*self.H6 + dWde*self.He + dWds*self.Hs
        return (
            self.mu * Identity2D,
            dWd4 * self.h4 + dWd6 * self.h6,
            dWde * self.He,
            dWds * self.Hs,
            array([self.mu]),
            array([self.c * dWd4 + self.c * dWd6]),
            I_n * Cinv,
        )
