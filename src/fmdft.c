/* Distributed 1D FFT using FMM

First each MPI process calls mfti to initialise parameters.
Then any number of calls can be made to mft to perform the transform using the same length and number of processors.
*/

#include <stdio.h>
#include <math.h>


/* Returns log to base 2 for an int */

int my_log2(int n) {
    int ans = 0;
    while(n>1) {
        n = n/2;
        ans++;
    }
    return ans;
}

void mftii() {

}

/*

mfti : Initialise the arrays for FMM-FFT

Inputs
`int lq` :  length of qr (total length / np)
`int np` :  number of procs
`int p` :  number of terms in expansions
`int nieach` :  number of particles in each box

Outputs
`double *w` :  work array of size at least:
`double *wt` :  work array of size at least:

where `n = log2(s)`, `lnp = log2(np)`, `s = lq/nieach`, `t = s/np`
*/
void mfti(int lq, int np, int p, int nieach, double *w, double *wt, int sz1, int szwk) {
    int s;
    /* ia is the index array */
    int ind, ia[30], j, n, lg2, nieach;

    s = lq / nieach;
    n = my_log2(s);
    ind = 1;

    ia[0] = ind;
    
    ind = ind + sz1;

    for(j=1; j<9; j++) {
        ia[j] = ind;
        ind += p*p*(n-2);
    }

    ia[9] = ind;
    ind += p*p;
    ia[10] = ind;
    ind += p*nieach;
    ia[11] = ind;
    ind += p*(np-1)*nieach;
    ia[12] = ind;
    ind += p*np/2;
    ia[13] = ind;
    ind += 3*nieach*np/2;
    ia[14] = ind;
    ind += 3*nieach*(nieach-1)*(np-1);
    ia[15] = ind;
    ind += 3*nieach*np/2;
    ia[16] = ind;
    ind += p;
    ia[17] = ind;
    ind += 2*(np-1)*np;
    ia[18] = ind;
    ind += nieach;
    ia[19] = ind;
    ind+=p;
    ia[20] = ind;
    ind+=p;
    ia[21] = ind;

    mftii(lq, np, s, p, nieach, n, sz1, szwk,
        w[ia[0]], w[ia[1]], w[ia[2]], w[ia[3]], w[ia[4]],
        w[ia[5]], w[ia[6]], w[ia[7]], w[ia[8]], w[ia[9]],
        w[ia[10]], w[ia[11]], w[ia[12]], w[ia[13]], w[ia[14]],
        w[ia[15]], w[ia[16]], w[ia[17]], w[ia[18]], w[ia[19]],
        w[ia[20]], w[ia[21]], wt
    );
}