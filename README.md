# Caputo Derivative

This is a companion source code for the publication Simulating Hyperelasticity and Fractional Viscoelasticity in the Human Heart.

## Content

This repository contains the following:
- Matlab caputo fractional derivative example (/matlab_example)
- Python caputo fractional derivative implementation (/src/py/caputoD)
- C++ caputo fractional derivative implementation (/src/cpp/fractional/caputo.cpp)
- Cython headers and module for compiling the c++ code to a python shared library (/src/cython)
- Script for generating the analytical solution in Example 1 of paper (compute_analytical_solutions.py)
- Script for generating the validation data for convergence with the number of Prony series (prony_serie_test.py)
- Script for generating the validation data for convergence with the scaling the time interval (time_interval_test.py)

## Installation

The necessary python modules are listed in requirements.txt
To compile the cython module, you can:
- Call Makefile (Linux only)
- Run python3 make.py
- Manually run pip install on requirements.txt and build /src/setup.py