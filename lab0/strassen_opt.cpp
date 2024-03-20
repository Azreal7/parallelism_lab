#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include <chrono>

using namespace std;

double* combineMatrix(double* C11, double* C12, double* C21, double* C22, int m, int n) {
    double* mat = (double*)malloc(m*n*8);
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            if(i < m / 2 && j < n / 2) mat[i * n + j] = C11[i * n / 2 + j];
            else if(i < m / 2 && j >= n / 2) mat[i * n + j] = C12[i * n / 2 + j - n / 2];
            else if(i >= m / 2 && j < n / 2) mat[i * n + j] = C21[(i - m / 2) * n / 2 + j];
            else mat[i * n + j] = C22[(i - m / 2) * n / 2 + j - n / 2];
        }
    }
    free(C11);
    free(C12);
    free(C21);
    free(C22);
    return mat;
}

double* ADD(double* mat1, double* mat2, int row, int col) {
    double* mat = (double*)malloc(row*col*8);
    for(int i = 0; i < row; i++) {
        for(int j = 0; j < col; j++) mat[i * col + j] = mat1[i * col + j] + mat2[i * col + j];
    }
    return mat;
}

double* SUB(double* mat1, double* mat2, int row, int col) {
    double* mat = (double*)malloc(row*col*8);
    for(int i = 0; i < row; i++) {
        for(int j = 0; j < col; j++) mat[i * col + j] = mat1[i * col + j] - mat2[i * col + j];
    }
    return mat;
}

double* MUL(double* mat1, double* mat2, int m, int n, int k) {
    if(m <= 128 || n <= 128 || k <= 128) {
        double* mat = (double*)malloc(m*k*8);
        memset(mat, 0, sizeof(mat));
        for(int i = 0; i < m; i++) {
            for(int l = 0; l < n; l++) {
                for(int j = 0; j < k; j++) {
                    mat[i * k + j] += mat1[i * n + l] * mat2[l * k + j];
                }
            }
        }
        return mat;
    }
    else {
        int m1 = m / 2, n1 = n / 2, k1 = k / 2;
        double* A11 = (double*)malloc(m1*n1*8);
        double* A12 = (double*)malloc(m1*n1*8);
        double* A21 = (double*)malloc(m1*n1*8);
        double* A22 = (double*)malloc(m1*n1*8);
        double* B11 = (double*)malloc(n1*k1*8);
        double* B12 = (double*)malloc(n1*k1*8);
        double* B21 = (double*)malloc(n1*k1*8);
        double* B22 = (double*)malloc(n1*k1*8);
        for(int i = 0; i < m1; i++) {
            for(int j = 0; j < n1; j++) {
                A11[i * n1 + j] = mat1[i * n + j];
                A12[i * n1 + j] = mat1[i * n + j + n1];
                A21[i * n1 + j] = mat1[(i + m1) * n + j];
                A22[i * n1 + j] = mat1[(i + m1) * n + j + n1];
            }
        }
        for(int i = 0; i < n1; i++) {
            for(int j = 0; j < k1; j++) {
                B11[i * k1 + j] = mat2[i * k + j];
                B12[i * k1 + j] = mat2[i * k + j + k1];
                B21[i * k1 + j] = mat2[(i + n1) * k + j];
                B22[i * k1 + j] = mat2[(i + n1) * k + j + k1];
            }
        }
        double* C11 = ADD(MUL(A11, B11, m1, n1, k1), MUL(A12, B21, m1, n1, k1), m1, k1);
        double* C12 = ADD(MUL(A11, B12, m1, n1, k1), MUL(A12, B22, m1, n1, k1), m1, k1);
        double* C21 = ADD(MUL(A21, B11, m1, n1, k1), MUL(A22, B21, m1, n1, k1), m1, k1);
        double* C22 = ADD(MUL(A21, B12, m1, n1, k1), MUL(A22, B22, m1, n1, k1), m1, k1);
        free(A11);
        free(A12);
        free(A21);
        free(A22);
        free(B11);
        free(B12);
        free(B21);
        free(B22);
        return combineMatrix(C11, C12, C21, C22, m, k);
    }
}

int main() {
    int m, n, k;
    scanf("%d%d%d", &m, &n, &k);
    double *mat1, *mat2;
    mat1 = (double*)malloc(m*n*8);
    mat2 = (double*)malloc(n*k*8);
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            mat1[i * n + j] = rand();
        }
    }
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < k; j++) {
            mat2[i * k + j] = rand();
        }
    }
    // printMatrix(mat1, m, n);
    auto start = chrono::high_resolution_clock::now();
    double *res = MUL(mat1, mat2, m, n, k);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "Original Time: " << duration.count() << " microseconds\n" << endl;
}