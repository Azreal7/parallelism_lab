#define MATSIZE 128
#include <random>
#include <stdio.h>
#include <algorithm>
#include <iostream>
#include "mpi.h"

using namespace std;

void matrix_multiply(int M, int N, int K, double *A, double *B, double *C) {
    for (int m = 0; m < M; m++) {
        for (int k = 0; k < K; k++) {
            C[m * K + k] = 0;
            for (int n = 0; n < N; n++) {
                C[m * K + k] += A[m * N + n] * B[n * K + k];
            }
        }
    }
}

void init_mat(double *mat, int row, int col) {
    for (int i = 0; i < row * col; ++i)
        mat[i] = rand();
}

int main() {
    // MPI初始化
    int comm_sz, my_rank;
    MPI_Init(nullptr, nullptr);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);
    MPI_Status status;
    int m = MATSIZE, n = MATSIZE, k = MATSIZE;

    int begin_row, end_row, avg_rows;
    avg_rows = k / comm_sz;
    double *mat1, *mat2, *mat, *local_mat1, *local_mat3;
    double start, end;
    mat2 = (double *)malloc(n * k * 8);
    if (my_rank == 0) {
        mat1 = (double *)malloc(m * n * 8);
        
        mat = (double *)malloc(m * k * 8);
        init_mat(mat1, m, n);
        init_mat(mat2, n, k);
        start = MPI_Wtime();
    }
    if (comm_sz == 1) {
        matrix_multiply(m, n, k, mat1, mat2, mat);
    } else {
        local_mat1 = (double *)malloc(avg_rows * n * 8);
        local_mat3 = (double *)malloc(avg_rows * k * 8);

        printf("#%d scattering.\n", my_rank);
        MPI_Scatter(mat1, avg_rows * n, MPI_DOUBLE, local_mat1, avg_rows * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        printf("#%d Bcasting.\n", my_rank);
        MPI_Bcast(mat2, n * k, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        matrix_multiply(avg_rows, n, k, local_mat1, mat2, local_mat3);

        printf("#%d Gathering.\n", my_rank);
        MPI_Gather(local_mat3, avg_rows * k, MPI_DOUBLE, mat, avg_rows * k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    if (my_rank == 0) {
        end = MPI_Wtime();
        printf("Total time: %lfs\n", (double)(end - start));
        free(mat1);
        free(mat);
    }
    if (comm_sz != 1) {
        free(local_mat1);
        free(local_mat3);
    }
    free(mat2);
    MPI_Finalize();
    return 0;
}