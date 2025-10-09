import dataclasses

import numpy as np

_FArray1D = np.ndarray[tuple[int], np.dtype[np.float64]]

Identity2D: _FArray1D = np.array([1.0, 0.0, 0.0, 1.0])


@dataclasses.dataclass(slots=True)
class Kinematics:
    det: float = 1.0
    I_n: float = 1.0
    I_1: float = 3.0
    I_1m3: float = 0.0
    C: _FArray1D = dataclasses.field(default_factory=lambda: np.array([1.0, 0.0, 0.0, 1.0]))
    Cinv: _FArray1D = dataclasses.field(default_factory=lambda: np.array([1.0, 0.0, 0.0, 1.0]))
    C33Cinv: _FArray1D = dataclasses.field(default_factory=lambda: np.array([1.0, 0.0, 0.0, 1.0]))


class DeformationGradient2D(Kinematics):
    def __init__(self, args: _FArray1D | None = None) -> None:
        super().__init__()
        if args is None:
            return
        self.precompute(args)

    def precompute(self, args: _FArray1D) -> None:
        self.det = args[0] * args[3] - args[1] * args[1]
        self.I_n = 1.0 / self.det
        self.C[:] = args[:]
        self.Cinv[0] = args[3] * self.I_n
        self.Cinv[1] = -args[1] * self.I_n
        self.Cinv[2] = -args[2] * self.I_n
        self.Cinv[3] = args[0] * self.I_n
        self.C33Cinv[:] = self.I_n * self.Cinv[:]
        self.I_1 = args[0] + args[3] + self.I_n
        self.I_1m3 = self.I_1 - 3.0
