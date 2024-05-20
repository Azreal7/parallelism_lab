#define MATSIZE 128
#include <pthread.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

int thread_num = 16;
int avg_size = MATSIZE / thread_num;

struct thread_data {
    int rank;
    double *A;
    double *B;
    double *C;
    thread_data(int _rank, double *_A, double *_B, double *_C) {
        rank = _rank;
        A = _A;
        B = _B;
        C = _C;
    }
};

// 串行化计算，用于与并行化结果进行比较
void is_equal(double *A, double *B, double *C) {
    for (int i = 0; i < MATSIZE; ++i) {
        for (int j = 0; j < MATSIZE; ++j) {
            double tmp = 0;
            for (int k = 0; k < MATSIZE; ++k) {
                tmp += A[i * MATSIZE + k] * B[k * MATSIZE + j];
            }
            if (abs(tmp - C[i * MATSIZE + j]) > 1e-6) {
                printf("false\n");
                return;
            }
        }
    }
    printf("true\n");
}

void mat_init(double *A, double *B) {
    for (int i = 0; i < MATSIZE; ++i) {
        for (int j = 0; j < MATSIZE; ++j) {
            A[i * MATSIZE + j] = random();
            B[i * MATSIZE + j] = random();
        }
    }
}

void print_mat(double *mat) {
    for (int i = 0; i < MATSIZE; ++i) {
        for (int j = 0; j < MATSIZE; ++j) {
            printf("%lf ", mat[i * MATSIZE + j]);
        }
        printf("\n");
    }
}

void *calculateUnit(void *arg) {
    int my_rank = ((thread_data *)arg)->rank;
    printf("It's thread %d\n", my_rank);
    int start_row = my_rank * avg_size;
    int end_row = (my_rank + 1) * avg_size;
    for (int i = start_row; i < end_row; ++i) {
        for (int j = 0; j < MATSIZE; ++j) {
            for (int k = 0; k < MATSIZE; ++k) {
                ((thread_data *)arg)->C[i * MATSIZE + j] +=
                    ((thread_data *)arg)->A[i * MATSIZE + k] * ((thread_data *)arg)->B[k * MATSIZE + j];
            }
        }
    }
    pthread_exit(nullptr);
}

int main() {
    double *A = (double *)malloc(MATSIZE * MATSIZE * sizeof(double));
    double *B = (double *)malloc(MATSIZE * MATSIZE * sizeof(double));
    double *C = (double *)malloc(MATSIZE * MATSIZE * sizeof(double));
    mat_init(A, B);
    thread_data *data = (thread_data *)malloc(thread_num * sizeof(thread_data));
    for (int i = 0; i < thread_num; ++i) {
        data[i] = thread_data(i, A, B, C);
    }
    pthread_t *threads;
    threads = (pthread_t *)malloc(thread_num * sizeof(pthread_t));

    // 开始时间记录
    auto start = std::chrono::steady_clock::now();

    for (int i = 0; i < thread_num; ++i) {
        pthread_create(&threads[i], nullptr, calculateUnit, (void *)&data[i]);
    }
    for (int i = 0; i < thread_num; ++i) {
        pthread_join(threads[i], nullptr);
    }

    // 结束时间记录
    auto end = std::chrono::steady_clock::now();

    // 计算时间
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    printf("Time: %lf s\n", (double)duration.count()/(double)1000000);
    // is_equal(A, B, C);
    return 0;
}