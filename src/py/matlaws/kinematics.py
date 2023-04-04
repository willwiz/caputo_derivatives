import numpy as np
import numpy.typing as npt
import typing as tp
import dataclasses

identity2D = np.array([1.0,0.0,0.0,1.0])

@dataclasses.dataclass
class kinematics:
    det:float = 1.0
    I_n:float = 1.0
    I_1:float = 3.0
    I_1m3:float = 0.0
    C:npt.NDArray[np.float64] = dataclasses.field(default_factory=lambda: np.array([1.0,0.0,0.0,1.0]))
    Cinv:npt.NDArray[np.float64] = dataclasses.field(default_factory=lambda: np.array([1.0,0.0,0.0,1.0]))
    C33Cinv:npt.NDArray[np.float64] = dataclasses.field(default_factory=lambda: np.array([1.0,0.0,0.0,1.0]))

class deformation_gradient_2D(kinematics):
    def __init__(self, args:tp.Optional[npt.NDArray[np.float64]]=None) -> None:
        super().__init__()
        if isinstance(args, np.ndarray):
            self.precompute(args)
    def precompute(self, args:tp.Optional[npt.NDArray[np.float64]]) -> None:
        self.det = args[0]*args[3] - args[1]*args[1]
        self.I_n = 1.0 / self.det
        self.C[:] = args[:]
        self.Cinv[0] =  args[3] * self.I_n
        self.Cinv[1] = -args[1] * self.I_n
        self.Cinv[2] = -args[2] * self.I_n
        self.Cinv[3] =  args[0] * self.I_n
        self.C33Cinv[:] = self.I_n * self.Cinv[:]
        self.I_1 = args[0] + args[3] + self.I_n
        self.I_1m3 = self.I_1 - 3.0