#include <pthread.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

int thread_num = 8;
pthread_mutex_t mtx;
double pi = 0;
int avg_num = 0;

void *thread_function(void *arg) {
    unsigned int seed = time(nullptr);
    int count = 0;
    for (int i = 0; i < avg_num; ++i) {
        double x = rand_r(&seed) / (double)RAND_MAX;
        double y = rand_r(&seed) / (double)RAND_MAX;
        if (x * x + y * y <= 1) {
            count += 1;
        }
    }
    pthread_mutex_lock(&mtx);
    pi += count;
    pthread_mutex_unlock(&mtx);
    pthread_exit(nullptr);
}

int my_main(int n) {
    avg_num = n / thread_num;
    pthread_t *threads = (pthread_t *)malloc(thread_num * sizeof(pthread_t));
    auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < thread_num; ++i) {
        pthread_create(&threads[i], nullptr, thread_function, nullptr);
    }
    for (int i = 0; i < thread_num; ++i) {
        pthread_join(threads[i], nullptr);
    }
    pi = 4 * pi / (double)n;
    printf("The estimated number of pi is %lf\n", pi);
    auto end = std::chrono::steady_clock::now();
    // 时间单位为毫秒
    printf("Time: %lf ms\n", std::chrono::duration<double, std::milli>(end - start).count());
    return 0;
}

int main() {
    int n[5] = {1024, 8192, 16364, 32768, 65536};
    for(int i = 0; i < 5; ++i) {
        my_main(n[i]);
    }
    return 0;
}