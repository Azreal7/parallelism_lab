#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include <chrono>
#define MINMATSIZE 128

using namespace std;

void ADD(double* mat1, double* mat2, double* mat, int col, int left, int right, int up, int down) {
    for(int i = left; i <= right; i++) {
        for(int j = up; j <= down; j++) mat[i * col + j] = mat1[i * col + j] + mat2[i * col + j];
    }
}

void SUB(double* mat1, double* mat2, double* mat, int col, int left, int right, int up, int down) {
    for(int i = left; i <= right; i++) {
        for(int j = up; j <= down; j++) mat[i * col + j] = mat1[i * col + j] - mat2[i * col + j];
    }
}

void MUL(double* mat1, double* mat2, double* mat, double *temp_mat, int m, int n, int k, int left1, int right1, int up1, int down1, int left2, int right2, int up2, int down2) {
    int m0 = right1 - left1, n0 = down1 - up1, k0 = down2 - up2;
    if(right1 - left1 <= 128 || up1 - down1 <= 128 || up2 - down2 <= 128) {
        for(int i = left1; i < right1; i++) {
            for(int q = up1; q < down1; q++) {
                for(int j = up2; j < down2; j++) {
                    mat[i * k + j] += mat1[i * n + q] * mat2[q * k + j];
                }
            }
        }
    }
    else {
        MUL(mat1, mat2, mat, temp_mat, m, n, k, left1, right1/2, up1, down1/2, left2, right2/2, up2, down2/2);
        MUL(mat1, mat2, temp_mat, mat, m, n, k, right1/2, right1, up1, down1/2, left2, right2/2, down2/2, down2);
        ADD()
    }
}

void test(double *mat1, double *mat2, double *mat, int m, int n, int k) {
    for(int i = 0; i < m; i++) {
        for(int l = 0; l < n; l++) {
            for(int j = 0; j < k; j++) {
                mat[i * k + j] += mat1[i * n + l] * mat2[l * k + j];
            }
        }
    }
}

bool cmp(double *mat1, double *mat2, int m, int k) {
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < k; j++) {
            if(mat1[i * k + j] != mat2[i * k + j]) return false;
        }
    }
    return true;
}

int main() {
    int m, n, k;
    scanf("%d%d%d", &m, &n, &k);
    double *mat1, *mat2, *mat, *temp_mat;
    mat1 = (double*)malloc(m*n*8);
    mat2 = (double*)malloc(n*k*8);
    mat = (double*)malloc(m*k*8);
    temp_mat = (double*)malloc(m*k*8);
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
    MUL(mat1, mat2, mat, m, n, k, 0, m, 0, n, 0, n, 0, k);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "Original Time: " << duration.count() << " microseconds\n" << endl;
    // test(mat1, mat2, temp_mat, m, n, k);
    // cout << cmp(temp_mat, mat, m, k) << endl;
}