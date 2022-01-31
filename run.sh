#!/bin/bash
rm output* time.txt
mpic++ -O3 amr.cpp -o amr
mpirun -np 4 amr > time.txt
