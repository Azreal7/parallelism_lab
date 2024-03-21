#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include <chrono>

using namespace std;

void printMatrix(double *mat, int row, int col) {
    for(int i = 0; i < row; i++) {
        for(int j = 0; j < col; j++) {
            printf("%.2lf ", mat[i * col + j]);
        }
        printf("\n");
    }
}

void multiplication(double *mat1, double *mat2, double *mat, int m, int n, int k) {
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < k; j++) {
            for(int l = 0; l < n; l++) {
                mat[i * k + j] += mat1[i * n + l] * mat2[l * k + j];
            }
        }
    }
}

void multiplication_loopUnwinding(double *mat1, double *mat2, double *mat, int m, int n, int k) {
    for(int i = 0; i < m; i++) {
        for(int l = 0; l < n; l++) {
            for(int j = 0; j < k; j+=8) {
                mat[i * k + j] += mat1[i * n + l] * mat2[l * k + j];
                mat[i * k + j+1] += mat1[i * n + l] * mat2[l * k + j+1];
                mat[i * k + j+2] += mat1[i * n + l] * mat2[l * k + j+2];
                mat[i * k + j+3] += mat1[i * n + l] * mat2[l * k + j+3];
                mat[i * k + j+4] += mat1[i * n + l] * mat2[l * k + j+4];
                mat[i * k + j+5] += mat1[i * n + l] * mat2[l * k + j+5];
                mat[i * k + j+6] += mat1[i * n + l] * mat2[l * k + j+6];
                mat[i * k + j+7] += mat1[i * n + l] * mat2[l * k + j+7];
            }
        }
    }
}

// 修改循环顺序
void new_multiplication(double *mat1, double *mat2, double *mat, int m, int n, int k) {
    for(int i = 0; i < m; i++) {
        for(int l = 0; l < n; l++) {
            for(int j = 0; j < k; j++) {
                mat[i * k + j] += mat1[i * n + l] * mat2[l * k + j];
            }
        }
    }
}

int main() {
    int m, n, k;
    // scanf("%d%d%d", &m, &n, &k);
    m = 1024, n = 1024, k = 1024;
    double *mat1, *mat2;
    mat1 = (double*)malloc(m*n*8);
    mat2 = (double*)malloc(n*k*8);
    double *mat = (double*)malloc(m*k*8);
    memset(mat, 0, sizeof(mat));
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
    auto start = chrono::high_resolution_clock::now();
    multiplication(mat1, mat2, mat, m, n, k);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "Original Time: " << duration.count() << " microseconds\n" << endl;

    auto start2 = chrono::high_resolution_clock::now();
    new_multiplication(mat1, mat2, mat, m, n, k);
    auto end2 = chrono::high_resolution_clock::now();
    auto duration2 = chrono::duration_cast<chrono::microseconds>(end2 - start2);
    cout << "New Time: " << duration2.count() << " microseconds\n" << endl;

    auto start1 = chrono::high_resolution_clock::now();
    multiplication_loopUnwinding(mat1, mat2, mat, m, n, k);
    auto end1 = chrono::high_resolution_clock::now();
    auto duration1 = chrono::duration_cast<chrono::microseconds>(end1 - start1);
    cout << "Loop unwinding Time: " << duration1.count() << " microseconds\n" << endl;
    // printMatrix(res, m, k);
}