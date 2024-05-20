#! /bin/bash
command g++ -fopenmp -o main main.cpp -L. -lparallelfor -lpthread
command export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.
command ./main
