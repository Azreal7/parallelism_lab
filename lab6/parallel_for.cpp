#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

struct thread_data {
    int start;
    int end;
    int ic;
    void *(*functor)(int, void*);
    void *arg;
};

void *thread_function(void *arg) {
    thread_data *data = (thread_data*)arg;
    int start = data->start;
    int end = data->end;
    int ic = data->ic;
    void *(*functor)(int, void*) = data->functor;
    void *t_arg = data->arg;
    for (int i = start; i < end; i += ic) {
        functor(i, t_arg);
    }
    pthread_exit(nullptr);
}


// is_circle表示是否是循环
void parallel_for(int start, int end, int ic, void *(*functor)(int, void*), void *arg, int num_threads, bool is_circle) {
    pthread_t threads[num_threads];
    struct thread_data datas[num_threads];

    // 计算线程任务范围
    int task_sum = (end - start + ic - 1) / ic;
    int avg_task = task_sum / num_threads; // 向下取整
    int r = task_sum % num_threads; // 余数

    for (int i = 0; i < num_threads; ++i) {
        if (!is_circle) {
            int task_start = start + i * avg_task * ic;
            int task_end = end;
            datas[i] = (struct thread_data){start + i * avg_task, start + (i + 1) * avg_task, ic, functor, arg};
        }
        else {
            int task_start = start + (i * avg_task + (i < r ? i : r)) * ic;
            int task_end = start + ((i + 1) * avg_task + (i <= r-1 ? i : r-1) + 1) * ic;
            datas[i] = (struct thread_data){task_start, task_end, ic, functor, arg};
        }
        pthread_create(&threads[i], nullptr, thread_function, &datas[i]);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], nullptr);
    }
}