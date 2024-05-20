#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

struct thread_data {
    int start;
    int end;
    int ic;
    void *(*functor)(int, void*);
    void *arg;
    thread_data(int start, int end, int ic, void *(*functor)(int, void*), void *arg) : start(start), end(end), ic(ic), functor(functor), arg(arg) {}
};

void *thread_function(void *arg) {
    thread_data *data = (thread_data*)arg;
    for (int i = data->start; i < data->end; i += data->ic) {
        data->functor(i, data->arg);
    }
    pthread_exit(nullptr);
}

void parallel_for(int start, int end, int ic, void *(*functor)(int, void*), void *arg, int num_threads) {
    pthread_t *threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    for (int i = 0; i < num_threads; ++i) {
        int avg_ic = (end - start) / num_threads;
        thread_data *data = (thread_data*)malloc(num_threads*sizeof(thread_data));
        data[i] = thread_data(start + i * avg_ic, start + (i + 1) * avg_ic, ic, functor, arg);
        pthread_create(&threads[i], nullptr, thread_function, (void*)data);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], nullptr);
    }
    free(threads);
}