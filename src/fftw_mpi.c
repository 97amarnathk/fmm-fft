/*******************************************************************************
*
*   Content : FFTW MPI C 1D complex, double precision
*
*******************************************************************************/

#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3-mpi.h>
#include <mpi.h>

int main(int argc, char* argv[]) {
    fftw_complex *x;

    int comm_size, comm_rank;
    int i, ret;
    fftw_plan forward_plan;
    ptrdiff_t N;
    ptrdiff_t size;
    ptrdiff_t local_ni, local_i_start, local_no, local_o_start;

    int RUNS = 1;
    double start, end;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    if(argc<2) {
        printf("Error give proper args!!");
        MPI_Finalize();
        return 0;
    }

    N = atoi(argv[1]);

    fftw_mpi_init();

    for(int run=0; run<RUNS; run++) {
        size = fftw_mpi_local_size_1d(N, MPI_COMM_WORLD, 
                FFTW_FORWARD, FFTW_ESTIMATE, &local_ni, 
                &local_i_start, &local_no, &local_o_start);

        x = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));

        for(int i=0; i<size; i++)
            x[i] = i;
        
        forward_plan = fftw_mpi_plan_dft_1d(N, x, x, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);

        start = MPI_Wtime();
        //--------------------------------------ALG STARTS HERE-----------------------------------
        fftw_execute(forward_plan);
        //--------------------------------------ALG ENDS  HERE-----------------------------------
        end = MPI_Wtime();

        if(0 == comm_rank) {
            printf("%td %d %d %lf\n", N, comm_size, run, end - start);
        }

        fftw_destroy_plan(forward_plan);
    }

    fftw_mpi_cleanup();
    MPI_Finalize();

    return 0;
}

/*
mpicc fftw_mpi.c -lm -lfftw3 -L/home/201501005/fftw3/fftw3/lib/ -I/home/201501005/fftw3/fftw3/include -lfftw3_mpi
gcc src/fftw_serial.c -lm -lfftw3 -L/home/201501005/fftw3/fftw3/lib/ -I/home/201501005/fftw3/fftw3/include -o fftw3_serial.o
*/