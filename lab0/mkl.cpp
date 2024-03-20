#include <mkl.h>
#include <iostream>
#include <chrono>

using namespace std;

int main() {
    double *A, *B, *C;
    int m, n, k, i, j;
    double alpha=1, beta=0;
    scanf("%d%d%d", &m, &n, &k);
    A = (double *)mkl_malloc( m*k*sizeof( double ), 64 );
    B = (double *)mkl_malloc( k*n*sizeof( double ), 64 );
    C = (double *)mkl_malloc( m*n*sizeof( double ), 64 );
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            A[i * n + j] = rand();
        }
    }
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < k; j++) {
            B[i * k + j] = rand();
        }
    }
    auto start = chrono::high_resolution_clock::now();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, k, B, n, beta, C, n);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "Time: " << duration.count() << " microseconds\n" << endl;
}
