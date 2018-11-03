#!/bin/bash

rm read
g++ geometry.C read.C analysis.C main.C -rpath ${ROOTSYS}/lib `root-config --libs --cflags` -o read # -rpath option might be necessary with some ROOT installations 
#compile for memory test
#g++ read.C `root-config --libs --cflags` -o read  -ltcmalloc
