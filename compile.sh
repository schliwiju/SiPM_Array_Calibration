#!/bin/bash

rm read
# -rpath option might be necessary with some ROOT installations 
# -lSpectrum option might be necessary with some ROOT installations 
g++ geometry.C read.C analysis.C main.C -rpath ${ROOTSYS}/lib `root-config --libs --cflags` -lSpectrum -o read 
#compile for memory test
#g++ read.C `root-config --libs --cflags` -o read  -ltcmalloc