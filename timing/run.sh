#!/bin/bash
#rm output* time.txt
#mpic++ -O3 amr.cpp -o amr
#mpirun -np 8 amr > time.txt
for i in 2 4 8 16 32
do
  filename="time"$i".txt"
  mpic++ -O3 amr.cpp -o amr
  mpirun -np $i amr > $filename
done
