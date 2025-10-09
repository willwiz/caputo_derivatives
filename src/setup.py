from pathlib import Path

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, setup

ext = [
    Extension(
        name="src.py.matlaws._h_models",
        sources=["src/cython/holzapfel_models.pyx"],
        include_dirs=[np.get_include(), str(Path(__file__).parents[1].resolve())],
    ),
]

setup(ext_modules=cythonize(ext, build_dir="src/cython/build"))
