#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <complex.h>

#define M_PI 3.14159265358979323846264338327950288

int my_log(int x) {
    int ans = 0;
    while(x>1) {
        x/=2;
        ans++;
    }
    return ans;
}

void initg(int b, int s, double* initch, int terms, double* wknie, double* chnods) {
    
}

void mftii(int lq, int p, int s, int terms, int b, int n, int szkeep, int sztemp,
        double* shftfn, double* shftfp, double* flip2n, double* flip2p,
        double* flip3n, double* flip3p, double* shftln, double* shftlp,
        double* topflp, double* initch, double* evalm, double* evalmh,
        double* cotsh, double* cots, double* cotprv, double chnods,
        double complex* fo, double* wknie, double* wkp, double* wkp3, double* wt)
    ) {
        double dlq = (double)lq;
        int t = s/p;
        b = lq/s;
        n = my_log2(s);

        for(int j=0; j<p; j++) {
            chnods[j] = cos((2*j-1)/(2*M_PI*terms));
        }

        initg(b, s, initch, terms, wknie, chnods);
        flip(2, 2, topflp, terms, chnods, wkp, wkp2, wkp3);

}

/*
initialise data for fmm-fft 
*/
void mfti(int lq, int p, int t, int b, double *wkkeep, double *wktemp, int szkeep, int sztemp) {
    // index array ia
    int ia[30], ind=0;
    int s = lq/b;
    int n = my_log2(s);

    ind = 0;

    for(int j=0; j<=7; j++) {
        ia[j] = ind;
        ind += t*t*(n-2);
    }

    ia[8] = ind;
    ind += t*t;
    ia[9] = ind;
    ind += t*b;
    ia[10] = ind;
    ind += t*(p-1)*b;
    ia[11] = ind;
    ind += (t*p)/2;
    ia[12] = ind;
    ind += (3*b*p)/2;
    ia[13] = ind;
    ind += 3*b*(b-1)*(p-1);
    ia[14] = ind;
    ind += (3*b*p)/2;
    ia[15] = ind;
    ind += t;    
    ia[16] = ind;
    ind+=2*(p-1)*p;
    ia[17] = ind;
    ind += b;

    for(int j=18; j<=20; j++) {
        ia[j] = ind;
        ind+=t;
    }

    /* Possible malloc on wkkeep required here */
    /*                                         */

    mftii(lq, p, s, t, b, n, szkeep, sztemp,
        &wkkeep[ia[0]], &wkkeep[ia[1]], &wkkeep[ia[2]], &wkkeep[ia[3]], &wkkeep[ia[4]],
        &wkkeep[ia[5]], &wkkeep[ia[6]], &wkkeep[ia[7]], &wkkeep[ia[8]], &wkkeep[ia[9]],
        &wkkeep[ia[10]], &wkkeep[ia[11]], &wkkeep[ia[12]], &wkkeep[ia[13]], &wkkeep[ia[14]],
        &wkkeep[ia[15]], &wkkeep[ia[16]], &wkkeep[ia[17]], &wkkeep[ia[18]], &wkkeep[ia[19]],
        &wkkeep[ia[20]], wktemp
    );
}

int main(int argc, char *argv[]) {
    //Initialise MPI env
    MPI_Init(&argc, &argv);
    int world_size, my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Initialise data using mfti
    mfti();

    //cleanup
    MPI_Finalize();
}