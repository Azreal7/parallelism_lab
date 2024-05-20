#! /bin/bash
command g++ -o matrix_mult matrix_mult.cpp -L. -lparallelfor -lpthread
command export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.
command ./matrix_mult
