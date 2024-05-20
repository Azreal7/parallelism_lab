#define MATSIZE 2048
#define NUM_THREADS 8
#include <omp.h>
#include <stdio.h>
#include <random>
#include <stdlib.h>
#include <chrono>

void init_mat(double *mat1, double *mat2) {
    for(int i = 0; i < MATSIZE; ++i) {
        for(int j = 0; j < MATSIZE; ++j) {
            mat1[i * MATSIZE + j] = (double)rand()/RAND_MAX;
            mat2[i * MATSIZE + j] = (double)rand()/RAND_MAX;
        }
    }
}

void is_equal(double *A, double *B, double *C) {
    for (int i = 0; i < MATSIZE; ++i) {
        for (int j = 0; j < MATSIZE; ++j) {
            double tmp = 0;
            for (int k = 0; k < MATSIZE; ++k) {
                tmp += A[i * MATSIZE + k] * B[k * MATSIZE + j];
            }
            if (abs(tmp - C[i * MATSIZE + j]) > 1e-6) {
                printf("%d %d %lf %lf\n", i, j, tmp, C[i * MATSIZE + j]);
                printf("false\n");
                return;
            }
        }
    }
    printf("true\n");
}

int main() {
    double *A = (double*)malloc(MATSIZE*MATSIZE*8);
    double *B = (double*)malloc(MATSIZE*MATSIZE*8);
    double *C = (double*)malloc(MATSIZE*MATSIZE*8);
    init_mat(A, B);
    int i, j, k;
    omp_set_num_threads(NUM_THREADS);
    auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel shared(A, B, C) private(j, k) 
    {
        #pragma omp for schedule(dynamic)
        for(i = 0; i < MATSIZE; ++i) {
            for(j = 0; j < MATSIZE; ++j) {
                C[i * MATSIZE + j] = 0;
                for(k = 0; k < MATSIZE; ++k) {
                    C[i * MATSIZE + j] += A[i * MATSIZE + k] * B[k * MATSIZE + j];
                }
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    printf("Time: %lf\n", elapsed.count());
    // is_equal(A, B, C);
}
