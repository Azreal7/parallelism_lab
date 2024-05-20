#define MATSIZE 1024
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

extern void parallel_for(int start, int end, int step, void*(*func)(int, void*), void* arg, int num_threads);

struct functor_args {
    double *A;
    double *B;
    double *C;
};

void init_mat(double *mat1, double *mat2) {
    for(int i = 0; i < MATSIZE; ++i) {
        for(int j = 0; j < MATSIZE; ++j) {
            mat1[i * MATSIZE + j] = (double)rand()/RAND_MAX;
            mat2[i * MATSIZE + j] = (double)rand()/RAND_MAX;
        }
    }
}

void *functor(int idx, void* args) {
    functor_args *fargs = (functor_args*)args;
    for (int j = 0; j < MATSIZE; ++j) {
        fargs->C[idx * MATSIZE + j] = 0;
        for (int k = 0; k < MATSIZE; ++k) {
            fargs->C[idx * MATSIZE + j] += fargs->A[idx * MATSIZE + k] * fargs->B[k * MATSIZE + j];
        }
    }
    return nullptr;
}

int my_main(int thread_num) {
    double *A = (double*)malloc(MATSIZE*MATSIZE*8);
    double *B = (double*)malloc(MATSIZE*MATSIZE*8);
    double *C = (double*)malloc(MATSIZE*MATSIZE*8);
    init_mat(A, B);
    functor_args args = {A, B, C};
    auto start = std::chrono::high_resolution_clock::now();
    parallel_for(0, MATSIZE, 1, functor, (void*)&args, thread_num);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    printf("Time: %lf\n", elapsed.count());
}

int main() {
    int nums[6] = {4,8,16,32,64,128};
    for (int i = 0; i < 6; ++i) {
        my_main(nums[i]);
    }
}