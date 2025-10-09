from pathlib import Path

import numpy as np


def write_array(file: Path | str, data: np.ndarray[tuple[int, ...], np.dtype[np.float64]]) -> None:
    with Path(file).open("w") as outfile:
        for i in data:
            outfile.writelines(f"{j:>22.12E}" for j in i)
            outfile.write("\n")


def write_list(file: Path | str, data: list[list[float]]) -> None:
    with Path(file).open("w") as outfile:
        for i in data:
            outfile.writelines(f"{j:>22}" for j in i)
            outfile.write("\n")


def write_array_1d(file: Path | str, data: np.ndarray[tuple[int], np.dtype[np.float64]]) -> None:
    with Path(file).open("w") as outfile:
        outfile.writelines(f"{i:>22.12E}" for i in data)


def read_array_2d_f(file: Path | str) -> np.ndarray[tuple[int, ...], np.dtype[np.float64]]:
    with Path(file).open() as outfile:
        arr = [[float(w) for w in line.strip().split()] for line in outfile if line.strip()]
    return np.array(arr)
