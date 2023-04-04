import numpy as np
import numpy.typing as npt
import typing as tp
from src.py.matlaws.kinematics import (
    identity2D,
    deformation_gradient_2D,
)

class NeoHookean_py:
    def __init__(self, mu:float) -> None:
        self.mu = mu
        pass
    def set_pars(self, mu:float) -> None:
        self.mu = mu
    def stress(self, strain:tp.Union[deformation_gradient_2D,npt.NDArray[np.float64]]):
        if isinstance(strain, deformation_gradient_2D):
            return self.mu * identity2D
        elif isinstance(strain, np.ndarray):
            kin = deformation_gradient_2D(strain)
            stress = self.stress(kin)
            return stress - self.mu * kin.C33Cinv
        else:
            raise TypeError('Arguments must by a vector of size 4 or deformation gradient class')