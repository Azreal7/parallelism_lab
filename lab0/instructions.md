mkl: g++ mkl.cpp -o mkl -m64  -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -m64  -I"${MKLROOT}/include" && ./mkl
original_c: g++ original_c.cpp -o original_c && ./original_c
original_c(编译优化): g++ original_c.cpp -o original_c -Ofast && ./original_c