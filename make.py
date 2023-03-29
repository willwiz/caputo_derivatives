#!/usr/bin/env python3

import os
import sys
import glob
import shutil
import subprocess

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
  empty_directory('src/py/*.pyd')
  empty_directory('src/py/*/*.pyd')
  empty_directory('src/py/*/*/*.pyd')
  empty_directory('src/py/*.so')
  empty_directory('src/py/*/*.so')
  empty_directory('src/py/*/*/*.so')


def main():
  print('Cleaning Build Folders')
  run_clean()
  print('\nChecking dependencies\n')
  subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-r', './requirements.txt'])
  print('\n')
  subprocess.check_call([sys.executable, './src/setup.py', 'build_ext', '--build-lib', '.'])