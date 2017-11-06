#!/bin/bash

cython cython_utils.pyx
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
      -I/usr/include/python2.7 -o cython_utils.so cython_utils.c

#python setup.py build_ext --inplace
