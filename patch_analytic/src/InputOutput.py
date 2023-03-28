#!/usr/bin/env python3
import numpy as np

def write_array(file, data):
    with open(file, 'w') as outfile:
        for i in data:
            for j in i:
                outfile.write('{:>22.12E}'.format(j))
            outfile.write('\n')
    return

def write_list(file, data):
    with open(file, 'w') as outfile:
        for i in data:
            for j in i:
                outfile.write('{:>22}'.format(j))
            outfile.write('\n')
    return

def write_array_1D(file, data):
    with open(file, 'w') as outfile:
        for i in data:
            outfile.write('{:>22.12E}'.format(i))
    return

def read_array_2D_f(file):
    with open(file, 'r') as outfile:
        arr = [[float(w) for w in line.strip().split()] for line in outfile if line.strip()]
    return np.array(arr)