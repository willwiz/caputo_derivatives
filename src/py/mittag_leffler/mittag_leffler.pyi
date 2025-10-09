import numpy as np

type _FArray1D = np.ndarray[tuple[int], np.dtype[np.float64]]

def ml[T: (_FArray1D, float, complex), F: (float, complex)](
    z: T, alpha: F, beta: F = 1.0, gamma: F = 1.0
) -> T: ...
def ml_[T: (_FArray1D, float, complex), F: (float, complex)](
    z: T, alpha: F, beta: F, gamma: F
) -> T: ...
