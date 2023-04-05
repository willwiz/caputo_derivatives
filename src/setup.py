from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import os
from pathlib import Path


ext = [
        Extension(name="src.py.matlaws._HModels",
          sources=["src/cython/HModels.pyx"],
          include_dirs=[numpy.get_include(), os.path.abspath(Path(__file__).parents[1])],
        ),
      ]

setup(ext_modules=cythonize(ext, build_dir="src/cython/build"))
