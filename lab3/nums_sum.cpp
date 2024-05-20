#include <pthread.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

int thread_num = 16;
int sum = 0;

struct thread_data {
    int rank;
    int *A;
    int n;
    thread_data(int _rank, int *_A, int _n) {
        rank = _rank;
        A = _A;
        n = _n;
    }
};

void *calculateUnit(void *arg) {
    int my_rank = ((thread_data *)arg)->rank;
    int *A = ((thread_data *)arg)->A;
    int avg_size = ((thread_data *)arg)->n / thread_num;
    int start = my_rank * avg_size;
    int end = (my_rank + 1) * avg_size;
    for (int i = start; i < end; ++i) {
        sum += A[i];
    }
    pthread_exit(NULL);
}

int main() {
    int ns[5] = {1, 8, 32, 64, 128};
    // scanf("%d", &n);
    for(int II = 0; II < 5; II++) {
        sum = 0;
        long long n = ns[II] * 1000000;
        int *A = (int*)malloc(n * sizeof(int));
        for(int i = 0; i < n; i++) {
            A[i] = rand();
        }
        pthread_t *threads = (pthread_t*)malloc(thread_num * sizeof(pthread_t));
        thread_data *data = (thread_data*)malloc(thread_num * sizeof(thread_data));
        for (int i = 0; i < thread_num; ++i) {
            data[i] = thread_data(i, A, n);
            pthread_create(&threads[i], NULL, calculateUnit, (void *)&data[i]);
        }
        auto start = std::chrono::steady_clock::now();
        for (int i = 0; i < thread_num; ++i) {
            pthread_join(threads[i], NULL);
        }
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        printf("Time: %lfs\n", elapsed_seconds.count());
    }
    
    return 0;
}