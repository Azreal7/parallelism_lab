#! /bin/bash
#
gcc -c -Wall -fopenmp heated_plate_openmp.c
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gcc -fopenmp heated_plate_openmp.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm heated_plate_openmp.o
mv a.out $HOME/parallelism_lab/lab6/heated_plate_openmp
#
echo "Normal end of execution."
