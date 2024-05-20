#define MATSIZE 128
#include <random>
#include <stdio.h>
#include <chrono>
#include <algorithm>
#include <iostream>
#include "mpi.h"

using namespace std;

template <typename T>
T fmax(T a, T b) {
    return a > b ? a : b;
}

void init_mat(double *mat, int row, int col) {
    for (int i = 0; i < row * col; i++)
        mat[i] = rand();
}

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

int main() {
    // MPI初始化
    int comm_sz, my_rank;
    MPI_Init(nullptr, nullptr);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);
    MPI_Status status;

    int m = MATSIZE, n = MATSIZE, k = MATSIZE;

    // 计算前n-2个进程计算的行的数量
    // 最后一个进程计算剩下大于avg_rows小于2*avg_rows行
    int begin_1row, end_1row, avg_rows;
    if (comm_sz > 1)
        avg_rows = k / (comm_sz - 1);
    if (my_rank == 0) {
        double *mat1, *mat2, *mat;
        mat1 = (double *)malloc(m * n * 8);
        mat2 = (double *)malloc(n * k * 8);
        mat = (double *)malloc(m * k * 8);
        init_mat(mat1, m, n);
        init_mat(mat2, n, k);
        auto start = chrono::high_resolution_clock::now();
        for (int i = 0; i < comm_sz - 1; ++i) {
            begin_1row = avg_rows * i,
            end_1row = (i + 1 == comm_sz - 1) ? fmax(k, avg_rows * (i + 1)) : avg_rows * (i + 1);
            MPI_Send(&end_1row, 1, MPI_INT, i + 1, 10, MPI_COMM_WORLD);
            MPI_Send(&mat1[begin_1row * n], (end_1row - begin_1row) * n, MPI_DOUBLE, i + 1, 1, MPI_COMM_WORLD);
            MPI_Send(mat2, n * k, MPI_DOUBLE, i + 1, 2, MPI_COMM_WORLD);
        }
        for (int i = 0; i < comm_sz - 1; ++i) {
            begin_1row = avg_rows * i,
            end_1row = (i + 1 == comm_sz - 1) ? fmax(k, avg_rows * (i + 1)) : avg_rows * (i + 1);
            MPI_Recv(&mat[begin_1row * n], (end_1row - begin_1row) * k, MPI_DOUBLE, i + 1, 3, MPI_COMM_WORLD, &status);
        }
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
        cout << "Time: " << duration.count() << " microseconds\n" << endl;
    } else {
        MPI_Recv(&end_1row, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &status);
        begin_1row = avg_rows * (my_rank - 1);

        printf("rank%d:\nfrom %d to %d\n", my_rank, begin_1row, end_1row);

        double *localMat1 = (double *)malloc((end_1row - begin_1row) * n * 8);
        double *localMat2 = (double *)malloc(n * k * 8);
        double *localMat3 = (double *)malloc((end_1row - begin_1row) * k * 8);

        MPI_Recv(localMat1, (end_1row - begin_1row) * n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(localMat2, n * k, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);

        matrix_multiply(end_1row - begin_1row, n, k, localMat1, localMat2, localMat3);

        MPI_Send(localMat3, (end_1row - begin_1row) * k, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);

        free(localMat1);
        free(localMat2);
        free(localMat3);
    }

    MPI_Finalize();
}