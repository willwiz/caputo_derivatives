from pathlib import Path
from typing import TypedDict

import numpy as np
from scipy.io import loadmat


class PrecomputedData(TypedDict):
    betam: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    taum: np.ndarray[tuple[int, int], np.dtype[np.float64]]


def import_precomputed_data(n_p: int) -> PrecomputedData:
    par_file = f"coeffs-opt-refined-100steps-{n_p}-500.mat"
    data = loadmat(Path(__file__).resolve().parent / "caputo_parm" / par_file)
    return {
        "betam": data["betam"],
        "taum": data["taum"],
    }
