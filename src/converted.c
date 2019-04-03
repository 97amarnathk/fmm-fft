#include <stdio.h>
#include <stdlib.h>
//#include <mpi.h>
#include <math.h>
#include <complex.h>

#define M_PI 3.14159265358979323846264338327950288

/* 
Returns 1D index given the 2D indices and array dimensions

Parameters:
`int i` : row index
`int j` : column index
`int dim1` : primary dimension of 2D array
`int dim2` : secondary dimension of 2D array

Returns:
1D index : (i*dim2 + j)
*/
int get1Dfrom2D(int i, int j, int dim1, int dim2) {
    return (i*dim2 + j);
}

/* 
Returns 1D index given the 3D indices and array dimensions

Parameters:
`int i` : primary index
`int j` : secondary index
`int k` : ternary index
`int dim1` : primary dimension of 3D array
`int dim2` : secondary dimension of 3D array
`int dim3` : ternary dimension of 3D array

Returns:
1D index : (i*dim2 + j)
*/
int get1Dfrom3D(int i, int j, int k, int dim1, int dim2, int dim3) {
    return((i * dim2 * dim3) + (j*dim3) + k);
}

int my_log(int x) {
    int ans = 0;
    while(x>1) {
        x/=2;
        ans++;
    }
    return ans;
}

void vx(double a, double* x, double* mat, int terms, double* chnods, double* mm, double* acu) {
    double ac = 3 * tan(a);

    for(int j=1; j<=terms; j++) {
        double th = -1 * cos((2*j-2)/(2 * (terms-1) * M_PI)) * (2*a);
        acu[j-1] = ac / tan(0.5 * th);
    }

    for(int k=1; k<=terms; k++) {
        mm[k-1] = chnods[k-1] - acu[k-1];

        for(int j=1; j<=k-1; j++) {
            mm[k-1] *= ((chnods[k-1] - acu[j-1])/(chnods[k-1]-chnods[j-1]));
        }

        for(int j=k+1; j<=terms; j++) {
            mm[k-1] *= ((chnods[k-1] - acu[j-1])/(chnods[k-1]-chnods[j-1]));
        }
    }

    for(int j=1; j<=terms; j++) {
        for(int k=1; k<=terms; k++) {
           mat[get1Dfrom2D(k-1, j-1, terms, terms)] = mm[j-1]/(x[k-1] - acu[j-1]);
           for(int l=1; l<=j-1; l++) {
               mat[get1Dfrom2D(k-1, j-1, terms, terms)] *= ((x[k-1] - chnods[l-1])/(x[k-1] - acu[l-1]));
           } 
           for(int l=j+1; l<=terms; l++) {
               mat[get1Dfrom2D(k-1, j-1, terms, terms)] *= ((x[k-1] - chnods[l-1])/(x[k-1] - acu[l-1]));  
           }
        }
    }

    return;
}

void flip(int lev, int shft, double* mat, int terms, double* chnods, double* x, double* wkp, double* mm) {
    double a = M_PI/(pow(2, lev));
    double b = shft;
    double c = tan(0.5 * a);
    double td2 = tan(b * a);

    for(int j=1; j<=p; j++) {
        x[j-1] = 3 * c * (1 - td2*c*chnods[j-1]) / (td2 + c*chnods[j-1]);
    }

    vx(0.5 * a, x, mat, terms, chnods, mm, wkp);

    return;
}

void initg(int b, int s, double* initch, int terms, double* tang, double* chnods) {
    double a = M_PI / (2.0 * s);

    for(int j=1; j<=b; j++) {
        th = ((double)(2*j - b - 1))/(double)b * a;
        tang[j-1] = tan(th);
    }

    double ta3 = 3 * tan(a);

    for(int j=1; j<=b; j++) {
        for(k=1; k<=p; k++) {
            initch[get1Dfrom2D(j-1, k-1, b, terms)] = (chnods[k-1] + ta3 * tang[j-1]) / (ta3 - chnods[k-1] * tang[j-1]);
        }
    }

    return;
}

void mftii(int lq, int p, int s, int terms, int b, int n, int szkeep, int sztemp,
        double* shftfn, double* shftfp, double* flip2n, double* flip2p,
        double* flip3n, double* flip3p, double* shftln, double* shftlp,
        double* topflp, double* initch, double* evalm, double* evalmh,
        double* cotsh, double* cots, double* cotprv, double* chnods,
        double complex* fo, double* wknie, double* wkp, double* wkp2, double* wkp3, double* wt
    ) {
        double dlq = (double)lq;
        int t = s/p;
        b = lq/s;
        n = my_log2(s);

        for(int j=1; j<=p; j++) {
            chnods[j-1] = cos((2*j-1)/(2*M_PI*terms));
        }

        initg(b, s, initch, terms, wknie, chnods);
        flip(2, 2, topflp, terms, chnods, wkp, wkp2, wkp3);

        for(int lev=1; lev<=n-2; lev++) {
            shftf(lev+2, -1, &shftfn[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2, wkp3);
            shftf(lev+2, 1, &shftfp[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2, wkp3);
            flip(lev+2, -2, &flip2n[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2, wkp3);
            flip(lev+2, +2, &flip2p[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2, wkp3);
            flip(lev+2, -3, &flip3n[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2, wkp3);
            flip(lev+2, +3, &flip3p[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2, wkp3);
            shftl(lev+1, -1, &shftln[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2);
            shftl(lev+1, 1, &shftlp[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2);
        }

        for(int sc=1; sc<=p-1; sc++) {
            double gfrac = ((double)sc)/p;

            for(int j=1; j<=b; j++) {
                double ffr = ((double)(2*j -1) - (((double)(2*sc))/p))/b - 1;
                evlmf(ffr, s, lq, evalm[get1Dfrom3D(j-1, sc-1, 0, b, p-1, terms)], terms, wkp, chnods, wkp2);
            }

            if(sc >= p/2) {
                int sce = sc - p/2 + 1;
                double ffr = (((double)2*b + 1) - ((double)2*sc)/(p))/(b) - 1;
                evlmf(ffr, s, lq, evalmh[get1Dfrom2D(sce-1, 0, p/2, terms)], terms, wkp, chnods, wkp2);

                for(int j= 1-b; j<=2*b; j++) {
                    cotsh[get1Dfrom2D(sce-1, j+b-1, p/2, 3*b)] = -1 / tan(M_PI/dlq * (j - 1 - b + gfrac)) / dlq;
                }
            }

            for(int k=2; k<=b; k++) {
                for(int j=1-b; j<=2*b; j++) {
                    cots[get1Dfrom3D(sc-1, k-2, j+b-1, p-1, b-1, 3*b)] = -1 / tan(M_PI/dlq * (j - k + gfrac)) / dlq;
                }
            }

            if(sc <= p/2) {
                for(int j=1-b; j<=2*b; j++) {
                    cotprv[get1Dfrom2D(sc-1, j + b -1, p/2, 3*b)] = -1 / tan(M_PI/dlq * (j - 1 + gfrac)) / dlq;
                }
            }
        }

        myffti();

        for(int sc = 1; sc<=p-1; sc++) {
            double ang = M_PI * (double)sc / p;
            for(int j=1; j<=p; j++) {
                //------> code to be added here :)
            }
        }
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