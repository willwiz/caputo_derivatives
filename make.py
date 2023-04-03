#!/usr/bin/env python3

import os
import sys
import glob
import shutil
import subprocess
import argparse


def remove_thing(path):
  if os.path.isdir(path):
    shutil.rmtree(path)
  else:
    os.remove(path)


def empty_directory(files):
  for i in glob.glob(files):
    remove_thing(i)


def run_clean(args=None):
  empty_directory('build/*')
  empty_directory('src/cython/build/*')
  empty_directory('src/py/**/*.pyd')
  empty_directory('src/py/**/*.so')


def run_build_dep():
  subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-r', './requirements.txt'])


def run_build():
  subprocess.check_call([sys.executable, './src/setup.py', 'build_ext', '--build-lib', '.'])


def run_all():
  run_clean()
  run_build_dep()
  run_build()


parser = argparse.ArgumentParser(description="""Makefile done with python""")
parser.set_defaults(main=run_build)
subparsers = parser.add_subparsers(help='sub-command help')
parser_all = subparsers.add_parser('all', help='a help')
parser_all.set_defaults(main=run_all)
parser_clean = subparsers.add_parser('clean', help='a help')
parser_clean.set_defaults(main=run_clean)
parser_build_dep = subparsers.add_parser('build-dep', help='a help')
parser_build_dep.set_defaults(main=run_build_dep)
parser_build = subparsers.add_parser('build', help='a help')
parser_build.set_defaults(main=run_build)


def main(args=None):
  args=parser.parse_args(args=args)
  args.main()


if __name__=="__main__":
  main()