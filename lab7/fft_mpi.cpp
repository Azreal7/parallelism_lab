#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <omp.h>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[]);
void ccopy(int n, double x[], double y[]);
void cfft2(int n, double x[], double y[], double w[], double sgn);
void cffti(int n, double w[]);
double cpu_time(void);
double ggl(double *ds);
void parallel_step(int n, int mj, double x[], double y[], double w[], double sgn);
void step(int n, int mj, double a[], double b[], double c[], double d[], double w[], double sgn);
void timestamp();

//****************************************************************************80

int main(int argc, char *argv[]) {
    int my_rank, comm_sz;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    double ctime;
    double ctime1;
    double ctime2;
    double error;
    int first;
    double flops;
    double fnm1;
    int i;
    int icase;
    int it;
    int ln2;
    double mflops;
    int n;
    int nits = 10000;
    static double seed;
    double sgn;
    double *w;
    double *x;
    double *y;
    double *z;
    double z0;
    double z1;

    if (my_rank == 0) {
        timestamp();
        cout << "\n";
        cout << "FFT_SERIAL\n";
        cout << "  C++ version\n";
        cout << "\n";
        cout << "  Demonstrate an implementation of the Fast Fourier Transform\n";
        cout << "  of a complex data vector.\n";
        //
        //  Prepare for tests.
        //
        cout << "\n";
        cout << "  Accuracy check:\n";
        cout << "\n";
        cout << "    FFT ( FFT ( X(1:N) ) ) == N * X(1:N)\n";
        cout << "\n";
        cout << "             N      NITS    Error         Time          Time/Call     MFLOPS\n";
        cout << "\n";
    }

    seed = 331.0;
    n = 1;
    //
    //  LN2 is the log base 2 of N.  Each increase of LN2 doubles N.
    //
    for (ln2 = 1; ln2 <= 20; ln2++) {
        n = 2 * n;
        //
        //  Allocate storage for the complex arrays W, X, Y, Z.
        //
        //  We handle the complex arithmetic,
        //  and store a complex number as a pair of doubles, a complex vector as a doubly
        //  dimensioned array whose second dimension is 2.
        //
        if (my_rank == 0) {
            w = new double[n];
            x = new double[2 * n];
            y = new double[2 * n];
            z = new double[2 * n];
        }

        first = 1;

        for (icase = 0; icase < 2; icase++) {
            if (my_rank == 0) {
                if (first) {
                    for (i = 0; i < 2 * n; i = i + 2) {
                        z0 = ggl(&seed);
                        z1 = ggl(&seed);
                        x[i] = z0;
                        z[i] = z0;
                        x[i + 1] = z1;
                        z[i + 1] = z1;
                    }
                } else {
                    for (i = 0; i < 2 * n; i = i + 2) {
                        z0 = 0.0;
                        z1 = 0.0;
                        x[i] = z0;
                        z[i] = z0;
                        x[i + 1] = z1;
                        z[i + 1] = z1;
                    }
                }
                //
                //  Initialize the sine and cosine tables.
                //
                cffti(n, w);
            }
            //
            //  Transform forward, back
            //
            if (first) {
                sgn = +1.0;
                cfft2(n, x, y, w, sgn);
                sgn = -1.0;
                cfft2(n, y, x, w, sgn);

                if (my_rank == 0) {
                    //
                    //  Results should be same as initial multiplied by N.
                    //
                    fnm1 = 1.0 / (double)n;
                    error = 0.0;
                    for (i = 0; i < 2 * n; i = i + 2) {
                        error = error + pow(z[i] - fnm1 * x[i], 2) + pow(z[i + 1] - fnm1 * x[i + 1], 2);
                    }
                    error = sqrt(fnm1 * error);
                    cout << "  " << setw(12) << n << "  " << setw(8) << nits << "  " << setw(12) << error;
                }
                first = 0;
            } else {
                ctime1 = cpu_time();
                for (it = 0; it < nits; it++) {
                    sgn = +1.0;
                    cfft2(n, x, y, w, sgn);
                    sgn = -1.0;
                    cfft2(n, y, x, w, sgn);
                }
                ctime2 = cpu_time();
                ctime = ctime2 - ctime1;

                flops = 2.0 * (double)nits * (5.0 * (double)n * (double)ln2);
                mflops = flops / 1.0E+06 / ctime;

                if (my_rank == 0) {
                    cout << "  " << setw(12) << ctime << "  " << setw(12) << ctime / (double)(2 * nits) << "  "
                         << setw(12) << mflops << "\n";
                }
            }
        }
        if ((ln2 % 4) == 0) {
            nits = nits / 10;
        }
        if (nits < 1) {
            nits = 1;
        }
        if (my_rank == 0) {
            delete[] w;
            delete[] x;
            delete[] y;
            delete[] z;
        }
    }

    if (my_rank == 0) {
        cout << "\n";
        cout << "FFT_SERIAL:\n";
        cout << "  Normal end of execution.\n";
        cout << "\n";
        timestamp();
    }

    MPI_Finalize();
    return 0;
}
//****************************************************************************80

void ccopy(int n, double x[], double y[]) {
    int i;

    for (i = 0; i < n; i++) {
        y[i * 2 + 0] = x[i * 2 + 0];
        y[i * 2 + 1] = x[i * 2 + 1];
    }
    return;
}
//****************************************************************************80

void cfft2(int n, double x[], double y[], double w[], double sgn) {
    int j;
    int m;
    int mj;
    int tgle;

    m = (int)(log((double)n) / log(1.99));
    mj = 1;
    //
    //  Toggling switch for work array.
    //
    tgle = 1;

    int my_rank, comm_sz;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    parallel_step(n, mj, x, y, w, sgn);

    if (n == 2) {
        return;
    }

    for (j = 0; j < m - 2; j++) {
        mj = mj * 2;
        if (tgle) {
            parallel_step(n, mj, y, x, w, sgn);
            tgle = 0;
        } else {
            parallel_step(n, mj, x, y, w, sgn);
            tgle = 1;
        }
    }
    //
    //  Last pass thru data: move y to x if needed
    //
    if (tgle) {
        if (my_rank == 0) {
            ccopy(n, y, x);
        }
    }

    mj = n / 2;
    parallel_step(n, mj, x, y, w, sgn);

    return;
}
//****************************************************************************80

void cffti(int n, double w[]) {
    double arg;
    double aw;
    int i;
    int n2;
    const double pi = 3.141592653589793;

    n2 = n / 2;
    aw = 2.0 * pi / ((double)n);

    for (i = 0; i < n2; i++) {
        arg = aw * ((double)i);
        w[i * 2 + 0] = cos(arg);
        w[i * 2 + 1] = sin(arg);
    }
    return;
}
//****************************************************************************80

double cpu_time(void) {
    double value;

    value = (double)clock() / (double)CLOCKS_PER_SEC;

    return value;
}
//****************************************************************************80

double ggl(double *seed) {
    double d2 = 0.2147483647e10;
    double t;
    double value;

    t = *seed;
    t = fmod(16807.0 * t, d2);
    *seed = t;
    value = (t - 1.0) / (d2 - 1.0);

    return value;
}
//****************************************************************************80

void parallel_step(int n, int mj, double x[], double y[], double w[], double sgn) {
    int my_rank, comm_sz;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Barrier(MPI_COMM_WORLD);
    if (comm_sz <= n / (2 * mj)) {
        int local_n = n / comm_sz;
        double *local_x = new double[2 * local_n];
        double *local_y = new double[2 * local_n];
        double *local_w = new double[local_n];

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Scatter(x, local_n, MPI_DOUBLE, local_x, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(x + n, local_n, MPI_DOUBLE, local_x + local_n, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(w, local_n, MPI_DOUBLE, local_w, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        step(local_n, mj, &local_x[0], &local_x[local_n], &local_y[0], &local_y[mj * 2], local_w, sgn);
        MPI_Gather(local_y, local_n * 2, MPI_DOUBLE, y, local_n * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        delete[] local_x;
        delete[] local_y;
        delete[] local_w;
    } else {  // n小于16的时候直接串行计算即可
        if (my_rank == 0) {
            step(n, mj, &x[0 * 2 + 0], &x[(n / 2) * 2 + 0], &y[0 * 2 + 0], &y[mj * 2 + 0], w, sgn);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
//****************************************************************************80

void step(int n, int mj, double a[], double b[], double c[], double d[], double w[], double sgn) {
    double ambr;
    double ambu;
    int j;
    int ja;
    int jb;
    int jc;
    int jd;
    int jw;
    int k;
    int lj;
    int mj2;
    double wjw[2];

    mj2 = 2 * mj;
    lj = n / mj2;

    for (j = 0; j < lj; j++) {
        jw = j * mj;
        ja = jw;
        jb = ja;
        jc = j * mj2;
        jd = jc;

        wjw[0] = w[(j * mj) * 2 + 0];
        wjw[1] = w[(j * mj) * 2 + 1];

        if (sgn < 0.0) {
            wjw[1] = -wjw[1];
        }

        for (k = 0; k < mj; k++) {
            c[(j * mj2 + k) * 2 + 0] = a[(j * mj + k) * 2 + 0] + b[(j * mj + k) * 2 + 0];
            c[(j * mj2 + k) * 2 + 1] = a[(j * mj + k) * 2 + 1] + b[(j * mj + k) * 2 + 1];

            ambr = a[(j * mj + k) * 2 + 0] - b[(j * mj + k) * 2 + 0];
            ambu = a[(j * mj + k) * 2 + 1] - b[(j * mj + k) * 2 + 1];

            d[(j * mj2 + k) * 2 + 0] = wjw[0] * ambr - wjw[1] * ambu;
            d[(j * mj2 + k) * 2 + 1] = wjw[1] * ambr + wjw[0] * ambu;
        }
    }
    return;
}
//****************************************************************************80

void timestamp() {
#define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    time_t now;

    now = time(NULL);
    tm = localtime(&now);

    strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

    cout << time_buffer << "\n";

    return;
#undef TIME_SIZE
}