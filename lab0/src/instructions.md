mkl: g++ mkl.cpp -o mkl -m64  -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -m64  -I"${MKLROOT}/include" && ./mkl
original_c: g++ original_c.cpp -o original_c -Ofast -mtune=native -march=native && ./original_c
strassen_opt: g++ strassen_opt.cpp -o strassen_opt -Ofast -mtune=native -march=native && ./strassen_opt