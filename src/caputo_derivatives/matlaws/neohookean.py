import numpy as np

from .kinematics import (
    DeformationGradient2D,
    Identity2D,
)

_FArray1D = np.ndarray[tuple[int], np.dtype[np.float64]]


class NeoHookeanPy:
    def __init__(self, mu: float) -> None:
        self.mu = mu

    def set_pars(self, mu: float) -> None:
        self.mu = mu

    def stress(self, strain: DeformationGradient2D | _FArray1D) -> _FArray1D:
        if isinstance(strain, DeformationGradient2D):
            return self.mu * Identity2D
        kin = DeformationGradient2D(strain)
        stress = self.stress(kin)
        return stress - self.mu * kin.C33Cinv
