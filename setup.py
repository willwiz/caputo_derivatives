# pyright: reportUnknownVariableType=false, reportUnknownArgumentType=false
from collections import defaultdict
from pathlib import Path

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

BUILD_ARGS: dict[str, list[str]] = defaultdict(
    lambda: ["-Ofast", "-g0", "-flto", "-march=native", "-std=c++17"]
)

BUILD_ARGS.update(
    {
        "msvc": ["/Ox"],
        "gcc": ["-Ofast", "-g0", "-flto", "-march=native", "-std=c++17"],
    }
)

include_dirs = [np.get_include(), "_cpp"]


class BuildExtCompilerCheck(build_ext):
    def build_extensions(self) -> None:
        compiler = self.compiler.compiler_type
        args = BUILD_ARGS[compiler]
        for ext in self.extensions:
            ext.extra_compile_args = args
            ext.include_dirs = include_dirs
        print(compiler)
        print(args)
        print(Path.cwd())
        build_ext.build_extensions(self)


ext = [
    Extension(
        name="caputo_derivatives.matlaws._h_models",
        sources=["_cython/holzapfel_models.pyx"],
    ),
]

print(Path().cwd())
setup(
    package_dir={"": "src"},
    ext_modules=cythonize(ext, build_dir="build"),
    cmdclass={"build_ext": BuildExtCompilerCheck},
)
