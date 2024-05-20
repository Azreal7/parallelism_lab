#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <pthread.h>

#define M           500
#define N           500
#define IS_CICLE    false

extern void parallel_for(
    int start,
    int end,
    int step,
    void *(*func)(int, void *),
    void *arg,
    int num_threads,
    bool is_cicle);

// 定义全局互斥锁
pthread_mutex_t mutex_for_mean = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_for_diff = PTHREAD_MUTEX_INITIALIZER;


// 定义functor参数结构体
struct functor_args_1 {
    double (*w)[N];
    int col;
    double value;
};

// 定义functor函数
// 作用：初始化W[idx][col] = value
void *functor_1(int idx, void *args) {
    struct functor_args_1 *args_data = (struct functor_args_1 *)args;
    args_data->w[idx][args_data->col] = args_data->value;
    return NULL;
}

// 定义functor参数结构体
struct functor_args_2 {
    double (*w)[N];
    int row;
    double value;
};

// 定义functor函数
// 作用：初始化W[row][idx] = value
void *functor_2(int idx, void *args) {
    struct functor_args_2 *args_data = (struct functor_args_2 *)args;
    args_data->w[args_data->row][idx] = args_data->value;
    return NULL;
}

// 定义functor参数结构体
struct functor_args_3 {
    double (*w)[N];
    double *mean;
};

// 定义functor函数
// 作用：mean = mean + w[idx][0] + w[idx][N-1];
void *functor_3(int idx, void *args) {
    struct functor_args_3 *args_data = (struct functor_args_3 *)args;

    pthread_mutex_lock(&mutex_for_mean);

    // 执行加法操作
    *(args_data->mean) = *(args_data->mean) + args_data->w[idx][0] + args_data->w[idx][N - 1];

    pthread_mutex_unlock(&mutex_for_mean);

    return NULL;
}

// 定义functor参数结构体
struct functor_args_4 {
    double (*w)[N];
    double *mean;
};

// 定义functor函数
// 作用：mean = mean + w[M-1][j] + w[0][j];
void *functor_4(int idx, void *args) {
    struct functor_args_4 *args_data = (struct functor_args_4 *)args;

    // 加锁
    pthread_mutex_lock(&mutex_for_mean);

    // 执行加法操作
    *(args_data->mean) = *(args_data->mean) + args_data->w[M - 1][idx] + args_data->w[0][idx];

    // 解锁
    pthread_mutex_unlock(&mutex_for_mean);

    return NULL;
}

// 定义functor参数结构体
struct functor_args_5 {
    double (*w)[N];
    double mean;
};

// 定义functor函数
// 作用：赋值w[idx][j] = mean
void *functor_5(int idx, void *args) {
    struct functor_args_5 *args_data = (struct functor_args_5 *)args;

    for (int j = 1; j < N - 1; j++) {
        args_data->w[idx][j] = args_data->mean;
    }

    return NULL;
}

// 定义functor参数结构体
struct functor_args_6 {
    double (*w)[N];
    double (*u)[N];
};

// 定义functor函数
// 作用：赋值u[idx][j] = w[idx][j]
void *functor_6(int idx, void *args) {
    struct functor_args_6 *args_data = (struct functor_args_6 *)args;

    for (int j = 0; j < N; j++) {
        args_data->u[idx][j] = args_data->w[idx][j];
    }

    return NULL;
}

// 定义functor参数结构体
struct functor_args_7 {
    double (*w)[N];
    double (*u)[N];
};

// 定义functor函数
// 作用：赋值w[idx][j] = (u[dix-1][j]+u[idx+1][j]+u[idx][j-1]+u[idx][j+1])/4.0
void *functor_7(int idx, void *args) {
    struct functor_args_7 *args_data = (struct functor_args_7 *)args;

    for (int j = 1; j < N - 1; j++) {
        args_data->w[idx][j] = (args_data->u[idx - 1][j] + args_data->u[idx + 1][j] + args_data->u[idx][j - 1] +
                                args_data->u[idx][j + 1]) /
            4.0;
    }

    return NULL;
}

// 定义functor参数结构体
struct functor_args_8 {
    double (*w)[N];
    double (*u)[N];
    double *my_diff;
};

// 定义functor函数
// 作用：my_diff = max(fabs(w[idx][j]-u[idx][j]))
void *functor_8(int idx, void *args) {
    struct functor_args_8 *args_data = (struct functor_args_8 *)args;

    for (int j = 1; j < N - 1; j++) {
        if (*(args_data->my_diff) < fabs(args_data->w[idx][j] - args_data->u[idx][j])) {
            *(args_data->my_diff) = fabs(args_data->w[idx][j] - args_data->u[idx][j]);
        }
    }

    return NULL;
}


int my_main(int num) {
    double diff;
    double epsilon = 0.001;
    int i;
    int iterations;
    int iterations_print;
    int j;
    double mean;
    double my_diff;
    double u[M][N];
    double w[M][N];
    double wtime;

    printf("\n");
    printf("HEATED_PLATE_OPENMP\n");
    printf("  C/OpenMP version\n");
    printf("  A program to solve for the steady state temperature distribution\n");
    printf("  over a rectangular plate.\n");
    printf("\n");
    printf("  Spatial grid of %d by %d points.\n", M, N);
    printf("  The iteration will be repeated until the change is <= %e\n", epsilon);
    printf("  Number of processors available = %d\n", omp_get_num_procs());
    printf("  Number of threads =              %d\n", omp_get_max_threads());

    /*
      Set the boundary values, which don't change.
    */
    mean = 0.0;

    struct functor_args_1 args_11 = {w, 0, 100.0};
    parallel_for(1, M - 1, 1, functor_1, (void *)&args_11, num, IS_CICLE);

    struct functor_args_1 args_12 = {w, N - 1, 100.0};
    parallel_for(1, M - 1, 1, functor_1, (void *)&args_12, num, IS_CICLE);

    struct functor_args_2 args_21 = {w, M - 1, 100.0};
    parallel_for(0, N, 1, functor_2, (void *)&args_21, num, IS_CICLE);

    struct functor_args_2 args_22 = {w, 0, 0.0};
    parallel_for(0, N, 1, functor_2, (void *)&args_22, num, IS_CICLE);

    /*
    Average the boundary values, to come up with a reasonable
    initial value for the interior.
    */
    struct functor_args_3 args_3 = {w, &mean};
    parallel_for(1, M - 1, 1, functor_3, (void *)&args_3, num, IS_CICLE);

    struct functor_args_4 args_4 = {w, &mean};
    parallel_for(0, N, 1, functor_4, (void *)&args_4, num, IS_CICLE);

    /*
      OpenMP note:
      You cannot normalize MEAN inside the parallel region.  It
      only gets its correct value once you leave the parallel region.
      So we interrupt the parallel region, set MEAN, and go back in.
    */
    mean = mean / (double)(2 * M + 2 * N - 4);
    printf("\n");
    printf("  MEAN = %f\n", mean);

    /*
      Initialize the interior solution to the mean value.
    */
    struct functor_args_5 args_5 = {w, mean};
    parallel_for(1, M - 1, 1, functor_5, (void *)&args_5, num, IS_CICLE);

    /*
      iterate until the  new solution W differs from the old solution U
      by no more than EPSILON.
    */
    iterations = 0;
    iterations_print = 1;
    printf("\n");
    printf(" Iteration  Change\n");
    printf("\n");
    wtime = omp_get_wtime();

    diff = epsilon;

    while (epsilon <= diff) {
        /*
          Save the old solution in U.
        */
        struct functor_args_6 args_6 = {w, u};
        parallel_for(0, M, 1, functor_6, (void *)&args_6, num, IS_CICLE);

        /*
          Determine the new estimate of the solution at the interior points.
          The new solution W is the average of north, south, east and west neighbors.
        */
        struct functor_args_7 args_7 = {w, u};
        parallel_for(1, M - 1, 1, functor_7, (void *)&args_7, num, IS_CICLE);

        /*
          C and C++ cannot compute a maximum as a reduction operation.

          Therefore, we define a private variable MY_DIFF for each thread.
          Once they have all computed their values, we use a CRITICAL section
          to update DIFF.
        */
        diff = 0.0;
        my_diff = 0.0;

        struct functor_args_8 args_8 = {w, u, &my_diff};
        parallel_for(1, M - 1, 1, functor_8, (void *)&args_8, num, IS_CICLE);

        pthread_mutex_lock(&mutex_for_diff);
        if (diff < my_diff) {
            diff = my_diff;
        }
        pthread_mutex_unlock(&mutex_for_diff);

        iterations++;
        if (iterations == iterations_print) {
            printf("  %8d  %f\n", iterations, diff);
            iterations_print = 2 * iterations_print;
        }
    }
    wtime = omp_get_wtime() - wtime;

    printf("\n");
    printf("  %8d  %f\n", iterations, diff);
    printf("\n");
    printf("  Error tolerance achieved.\n");
    printf("  Wallclock time = %f\n", wtime);
    /*
      Terminate.
    */
    printf("\n");
    printf("HEATED_PLATE_OPENMP:\n");
    printf("  Normal end of execution.\n");

    return 0;

#undef M
#undef N
}

int main() {
    int nums[5] = {1, 2, 4, 8, 16};
    for (int i = 0; i < 5; i++) {
        my_main(nums[i]);
        printf("%d threads finished\n", nums[i]);
    }
    return 0;
}