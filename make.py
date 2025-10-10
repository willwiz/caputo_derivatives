import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from warnings import warn


def remove_thing(path: Path | str) -> None:
    path = Path(path)
    if path.is_dir():
        shutil.rmtree(path)
        return
    if path.is_file():
        path.unlink()
    warn(f"Could not find {path} to remove", RuntimeWarning, stacklevel=2)


def empty_directory(files: Path | str) -> None:
    files = Path(files)
    for i in files.glob("*"):
        remove_thing(i)


def run_clean() -> None:
    empty_directory("build/*")
    empty_directory("src/cython/build/*")
    empty_directory("src/**/*.pyd")
    empty_directory("src/**/*.so")


def run_build_dep() -> None:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", "./requirements.txt"])


def run_build() -> None:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "."])


def run_test() -> None:
    subprocess.check_call([sys.executable, "-m", "unittest", "discover", "./tests"])


def run_all() -> None:
    run_clean()
    run_build_dep()
    run_build()
    run_test()


parser = argparse.ArgumentParser(description="""Makefile done with python""")
parser.set_defaults(main=run_build)
subparsers = parser.add_subparsers(help="sub-command help")
parser_all = subparsers.add_parser("all", help="run clean build-dep build test")
parser_all.set_defaults(main=run_all)
parser_clean = subparsers.add_parser("clean", help="clean shared lib and build files")
parser_clean.set_defaults(main=run_clean)
parser_build_dep = subparsers.add_parser("build-dep", help="install python dependencies from pip")
parser_build_dep.set_defaults(main=run_build_dep)
parser_build = subparsers.add_parser("build", help="compile cython shared lib")
parser_build.set_defaults(main=run_build)
parser_test = subparsers.add_parser("test", help="run unit tests")
parser_test.set_defaults(main=run_test)


def main(c_args: list[str] | None = None) -> None:
    args = parser.parse_args(args=c_args)
    args.main()


if __name__ == "__main__":
    main()
