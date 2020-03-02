#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <mpi.h>
#include <omp.h>

#include "fftw3.h"

/* Compute (K*L)%M accurately */
static double moda(int K, int L, int M)
{
    return (double)(((long long)K * L) % M);
}

/* Initialize array x[N] with harmonic H */
static void init(fftw_complex *x, int N, int H)
{
    double TWOPI = 6.2831853071795864769, phase;
    int n;

    for (n = 0; n < N; n++)
    {
        phase  = moda(n,H,N) / N;
        x[n][0] = cos( TWOPI * phase ) / N;
        x[n][1] = sin( TWOPI * phase ) / N;
    }
}

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    if(argc < 2) {
        printf("Error : give args\n");
        return 0;
    }

    int RUNS=6;

    int N = atoi(argv[1]);
    int i=0;

    int H =  -N/2;
    fftw_plan forward_plan = 0, backward_plan = 0;
    fftw_complex *x = 0;
    int status = 0;

    x = fftw_malloc(sizeof(fftw_complex)*N);

    for(i=0; i<RUNS; i++) {
        int runId = i;
        forward_plan = fftw_plan_dft(1, &N, x, x, FFTW_FORWARD, FFTW_ESTIMATE);

        // double start_time = MPI_Wtime();
        double start_time = omp_get_wtime();
        /*--------------ALG STARTS HERE --------------------------*/
        fftw_execute(forward_plan);
        /*--------------ALG ENDS HERE --------------------------*/
        double end_time = omp_get_wtime() - start_time;
        // double end_time = MPI_Wtime();

        printf("%d, %d, %lf\n", N, runId, end_time);

        // fftw_cleanup_threads();
        fftw_destroy_plan(forward_plan);
    }

    fftw_free(x);
    MPI_Finalize();
}
