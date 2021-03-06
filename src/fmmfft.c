/*******************************************************************************
*
*   Content : FMM-FFT MPI C 1D complex, double precision.
*   Authors :
*       Amarnath Karthi  (201501005@daiict.ac.in)
*       Chahak Mehta     (201501422@daiict.ac.in)
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <complex.h>
#include <omp.h>

#include <fftw3.h>

double PI  = 3.14159265358979323846264338327950288;

/* 
Returns 1D index given the 2D indices and array dimensions

Parameters:
`long int i` : row index
`long int j` : column index
`long int dim1` : primary dimension of 2D array
`long int dim2` : secondary dimension of 2D array

Returns:
1D index : (i*dim2 + j)
*/
long int get1Dfrom2D(long int i, long int j, long int dim1, long int dim2) {
    return (i*dim2 + j);
}

/* 
Returns 1D index given the 3D indices and array dimensions

Parameters:
`long int i` : primary index
`long int j` : secondary index
`long int k` : ternary index
`long int dim1` : primary dimension of 3D array
`long int dim2` : secondary dimension of 3D array
`long int dim3` : ternary dimension of 3D array

Returns:
1D index : (i*dim2 + j)
*/
long int get1Dfrom3D(long int i, long int j, long int k, long int dim1, long int dim2, long int dim3) {
    return((i * dim2 * dim3) + (j*dim3) + k);
}

long int get1Dfrom4D(long int i, long int j, long int k, long int l, long int dim1, long int dim2, long int dim3, long int dim4) {
    return((i * dim2 * dim3 * dim4) + (j*dim3*dim4) + k*dim4 + l);
}

long int my_log2(long int x) {
    long int ans = 0;
    while(x>1) {
        x/=2;
        ans++;
    }
    return ans;
}

void cjall(long int n, double complex* x) {
    for(long int i=0; i<n; i++) {
        x[i] = conj(x[i]);
    }
}

double complex dotp(long int n, double complex* x, double* r) {
    double complex ans = 0;
    for(long int i=0; i<n; i++) {
        ans += x[i]*r[i];
    }

    return ans;
}

void vx(double a, double* x, double* mat, long int terms, double* chnods, double* mm, double* acu) {
    double ac = 3 * tan(a);

    for(long int j=1; j<=terms; j++) {
        double th = -1 * cos(((double)2*j-2)/(2 * (terms-1)) * PI) * (2*a);
        acu[j-1] = ac / tan(0.5 * th);
    }

    for(long int k=1; k<=terms; k++) {
        mm[k-1] = chnods[k-1] - acu[k-1];

        for(long int j=1; j<=k-1; j++) {
            mm[k-1] *= ((chnods[k-1] - acu[j-1])/(chnods[k-1]-chnods[j-1]));
        }

        for(long int j=k+1; j<=terms; j++) {
            mm[k-1] *= ((chnods[k-1] - acu[j-1])/(chnods[k-1]-chnods[j-1]));
        }
    }

    for(long int j=1; j<=terms; j++) {
        for(long int k=1; k<=terms; k++) {
           mat[get1Dfrom2D(k-1, j-1, terms, terms)] = mm[j-1]/(x[k-1] - acu[j-1]);
           for(long int l=1; l<=j-1; l++) {
               mat[get1Dfrom2D(k-1, j-1, terms, terms)] *= ((x[k-1] - chnods[l-1])/(x[k-1] - acu[l-1]));
           } 
           for(long int l=j+1; l<=terms; l++) {
               mat[get1Dfrom2D(k-1, j-1, terms, terms)] *= ((x[k-1] - chnods[l-1])/(x[k-1] - acu[l-1]));  
           }
        }
    }

    return;
}

void wx(double a, double* x, long int xlen, double* mat, long int terms, double* chnods, double* acu) {
    double ac = 1/tan(a);

    for(long int j=1; j<=terms; j++) {
        double th = -cos( ((double)(2*j - 2)) / (2*(terms-1)) * PI ) * (PI - 6*a) + PI;
        acu[j-1] = ac * tan(0.5 * th);
        chnods[j-1] = -cos(((double)2*j -1) / (2*terms) * PI);
    }

    for(long int j=1; j<=terms; j++) {
        double mm = chnods[j-1] - acu[j-1];

        for(long int k=1; k<=j-1; k++) {
            mm *= ( (chnods[j-1] - acu[k-1]) / (chnods[j-1] - chnods[k-1]) );
        }

        for(long int k=j+1; k<=terms; k++) {
            mm *= ( (chnods[j-1] - acu[k-1]) / (chnods[j-1] - chnods[k-1]) );
        }
        
        for(long int k=1; k<=xlen; k++) {
            mat[get1Dfrom2D(k-1, j-1, xlen, terms)] = mm / (x[k-1] - acu[j-1]);
            for(long int l=1; l<=j-1; l++) {
                mat[get1Dfrom2D(k-1, j-1, xlen, terms)] *= ( (x[k-1] - chnods[l-1]) / (x[k-1] - acu[l-1]) );
            }
            for(long int l=j+1; l<=terms; l++) {
                mat[get1Dfrom2D(k-1, j-1, xlen, terms)] *= ( (x[k-1] - chnods[l-1]) / (x[k-1] - acu[l-1]) );
            }
        }
    }
}

void shftf(long int lev, long int dir, double* mat, long int terms, double* chnods, double* x, double* wkp, double* mm) {
    double dd = (double) dir;
    double a = PI/pow(2, lev);
    double ta3 = 3 * tan(a);
    double td2 = tan(0.5 * a);

    for(long int j=1; j<=terms; j++) {
        x[j-1] = 3 * td2 * (chnods[j-1] + dd*ta3*td2) / (ta3 - dd*td2*chnods[j-1]);
    }

    vx(0.5*a, x, mat, terms, chnods, wkp, mm);
}

void evlmf(double ffr, long int s, long int lq, double* vec, long int terms, double* w, double* chnods, double* acu) {
    double a = (PI / s) * 0.5;
    double th = ffr * a;
    double x = tan(th)/tan(a);
    wx(a, &x, 1, w, terms, chnods, acu);

    for(long int t = 1; t<=terms; t++) {
        vec[t-1] = w[t-1]/lq;
    }

    return;
}

void mpp(double complex* vecs, long int terms, long int len, double* matpp, double complex* ans) {
    for(long int sc=1; sc<=len; sc++) {
        for(long int term=1; term<=terms; term++) {
            ans[get1Dfrom2D(sc-1, term-1,len, terms)] = vecs[get1Dfrom2D(sc-1, 0,len ,terms)] * matpp[get1Dfrom2D(term-1, 0, terms, terms)];

            for(long int k=2; k<=terms; k++) {
                ans[get1Dfrom2D(sc-1, term-1,len, terms)] += vecs[get1Dfrom2D(sc-1, k-1,len ,terms)] * matpp[get1Dfrom2D(term-1, k-1, terms, terms)];
            }
        }
    }
}

void shftl(long int lev, long int dir, double* mat, long int terms, double* chnods, double* x, double* acu) {
    double dd = (double)dir;
    double a = PI/pow(2, lev);
    double c = tan(0.25 * a);
    double td2 = tan(0.5 * a);

    for(long int j=1; j<=terms; j++) {
        x[j-1] = (chnods[j-1] - dd) / (1/c + dd*c*chnods[j-1]) / td2;
    }

    wx(0.5 * a, x, terms, mat, terms, chnods, acu);
}

void initmm(double complex* qr, long int p, long int b, long int t, double* initch, double complex* phi, long int terms) {
    for(long int box=1; box<=t; box++) {
        for(long int sc=1; sc<=p-1; sc++) {
            for(long int term=1; term<=terms; term++) {
                phi[get1Dfrom3D(box-1, sc-1, term-1, t, p-1, terms)] = qr[get1Dfrom3D(box-1, sc, 0, t, p, b)] * initch[get1Dfrom2D(0, term-1, b, terms)];

                for(long int j=2; j<=b; j++) {
                    phi[get1Dfrom3D(box-1, sc-1, term-1, t, p-1, terms)] += qr[get1Dfrom3D(box-1, sc, j-1, t, p, b)] * initch[get1Dfrom2D(j-1, term-1, b, terms)];
                }
            }
        }
    }
}

double complex mn3(long int n, double complex* x, double complex* y, double complex* z, double* r) {
    double complex ans = x[0] * r[0];
    for(long int j=2; j<=n; j++) {
        ans += x[j-1] * r[j-1];
    }
    for(long int j=1; j<=n; j++) {
        ans+= y[j-1]*r[n+j-1];
    }
    for(long int j=1; j<=n; j++) {
        ans+= z[j-1]*r[n+n+j-1];
    }
    return ans;
}

void flip(long int lev, long int shft, double* mat, long int terms, double* chnods, double* x, double* wkp, double* mm) {
    double a = PI/(pow(2, lev));
    double b = shft;
    double c = tan(0.5 * a);
    double td2 = tan(b * a);

    for(long int j=1; j<=terms; j++) {
        x[j-1] = 3 * c * (1 - td2*c*chnods[j-1]) / (td2 + c*chnods[j-1]);
    }

    vx(0.5 * a, x, mat, terms, chnods, mm, wkp);

    return;
}

void initg(long int b, long int s, double* initch, long int terms, double* tang, double* chnods) {
    double a = PI / (2.0 * s);

    for(long int j=1; j<=b; j++) {
        double th = ((double)(2*j - b - 1))/(double)b * a;
        tang[j-1] = tan(th);
    }

    double ta3 = 3 * tan(a);

    for(long int j=1; j<=b; j++) {
        for(long int k=1; k<=terms; k++) {
            initch[get1Dfrom2D(j-1, k-1, b, terms)] = (chnods[k-1] + ta3 * tang[j-1]) / (ta3 - chnods[k-1] * tang[j-1]);
        }
    }

    return;
}

void mftii(long int lq, long int p, long int s, long int terms, long int b, long int n, long int szkeep, long int sztemp,
        double* shftfn, double* shftfp, double* flip2n, double* flip2p,
        double* flip3n, double* flip3p, double* shftln, double* shftlp,
        double* topflp, double* initch, double* evalm, double* evalmh,
        double* cotsh, double* cots, double* cotprv, double* chnods,
        double complex* fo, double* wknie, double* wkp, double* wkp2, double* wkp3, double* wt
    ) {
        double dlq = (double)lq;
        long int t = s/p;
        b = lq/s;
        n = my_log2(s);

        for(long int j=1; j<=terms; j++) {
            chnods[j-1] = -1*cos(((double)2*j-1)/(2*terms) * PI);
        }

        initg(b, s, initch, terms, wknie, chnods);

        flip(2, 2, topflp, terms, chnods, wkp, wkp2, wkp3);

        for(long int lev=1; lev<=n-2; lev++) {
            shftf(lev+2, -1, &shftfn[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2, wkp3);
            shftf(lev+2, 1, &shftfp[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2, wkp3);
            flip(lev+2, -2, &flip2n[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2, wkp3);
            flip(lev+2, +2, &flip2p[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2, wkp3);
            flip(lev+2, -3, &flip3n[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2, wkp3);
            flip(lev+2, +3, &flip3p[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2, wkp3);
            shftl(lev+1, -1, &shftln[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2);
            shftl(lev+1, 1, &shftlp[get1Dfrom3D(lev-1, 0, 0, n-2, terms, terms)], terms, chnods, wkp, wkp2);
        }

        for(long int sc=1; sc<=p-1; sc++) {
            double gfrac = ((double)sc)/p;

            for(long int j=1; j<=b; j++) {
                double ffr = ((double)(2*j -1) - (((double)(2*sc))/p)) / b - 1;
                evlmf(ffr, s, lq, &evalm[get1Dfrom3D(j-1, sc-1, 0, b, p-1, terms)], terms, wkp, chnods, wkp2);
            }

            if(sc >= p/2) {
                long int sce = sc - p/2 + 1;
                double ffr = (((double)2*b + 1) - ((double)2*sc)/(p))/(b) - 1;
                evlmf(ffr, s, lq, &evalmh[get1Dfrom2D(sce-1, 0, p/2, terms)], terms, wkp, chnods, wkp2);

                for(long int j= 1-b; j<=2*b; j++) {
                    cotsh[get1Dfrom2D(sce-1, j+b-1, p/2, 3*b)] = -1 / tan(PI/dlq * (j - 1 - b + gfrac)) / dlq;
                }
            }

            for(long int k=2; k<=b; k++) {
                for(long int j=1-b; j<=2*b; j++) {
                    cots[get1Dfrom3D(sc-1, k-2, j+b-1, p-1, b-1, 3*b)] = -1 / tan(PI/dlq * (j - k + gfrac)) / dlq;
                }
            }

            if(sc <= p/2) {
                for(long int j=1-b; j<=2*b; j++) {
                    cotprv[get1Dfrom2D(sc-1, j + b -1, p/2, 3*b)] = -1 / tan(PI/dlq * (j - 1 + gfrac)) / dlq;
                }
            }
        }

        //myffti();

        for(long int sc = 1; sc<=p-1; sc++) {
            double ang = PI * ((double)sc) / p;
            for(long int j=1; j<=p; j++) {
                double complex temp = (ang * (1 - 2*j))*I;
                fo[get1Dfrom2D(j-1, sc-1, p, p-1)] = -cexp(temp)*sin(ang);
            }
        }
}

/*
initialise data for fmm-fft 
*/ 
void mfti(long int lq, long int p, long int t, long int b, double *wkkeep, double *wktemp, long int szkeep, long int sztemp) {
    // index array ia
    long int ia[30], ind=0;
    long int s = lq/b;
    long int n = my_log2(s);

    ind = 0;

    for(long int j=0; j<=7; j++) {
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

    for(long int j=18; j<=20; j++) {
        ia[j] = ind;
        ind+=t;
    }
    mftii(lq, p, s, t, b, n, szkeep, sztemp,
        &wkkeep[ia[0]], &wkkeep[ia[1]], &wkkeep[ia[2]], &wkkeep[ia[3]], &wkkeep[ia[4]],
        &wkkeep[ia[5]], &wkkeep[ia[6]], &wkkeep[ia[7]], &wkkeep[ia[8]], &wkkeep[ia[9]],
        &wkkeep[ia[10]], &wkkeep[ia[11]], &wkkeep[ia[12]], &wkkeep[ia[13]], &wkkeep[ia[14]],
        &wkkeep[ia[15]], &wkkeep[ia[16]], &wkkeep[ia[17]], &wkkeep[ia[18]], &wkkeep[ia[19]],
        &wkkeep[ia[20]], wktemp
    );
}

void mftint(double complex* qr, long int lq, long int p, long int myid, long int s, long int terms, long int n, long int t, long int b, long int log2np, long int sz1, long int szwk, long int dir,
    double* shftfn, double* shftfp, double* flip2n, double* flip2p, double* flip3n, double* flip3p,
    double* shftln, double* shftlp, double* topflp, double* initch, double* evalm, double* evalmh,
    double* cotsh, double* cots, double* cotprv, double* cnods, double complex* fo,
    double complex* phi, double complex* psi, double complex* sump, double complex* sumsec, double complex* phils, double complex* phirs, double complex* phr,
    double complex* packps, double complex* packnr, double complex* packns, double complex* packpr, double complex* qrn, double complex* qrp, double complex* qcpp, double complex* packnn,
    double complex* phiopp, double complex* phim2, double complex* phip2,double complex* phim3, double complex* phip3, double complex* shpsi, double complex* psiprv,
    double complex* f2n, double complex* f2p, double complex* f3, double complex* phin, double complex* phip, double complex* psiev, double complex* psit, double complex* exts, double complex* extr,
    double* wknie, double* wkp, double* wkp2, double* wkp3,
    long int* prevd, long int* nextd
) {

    MPI_Barrier(MPI_COMM_WORLD);

    if(dir<0) {
        cjall(lq, qr);
    }

    double dlq = (double)lq;
    log2np = my_log2(p);
    long int npm1 = p-1;
    long int pnpm1 = terms*npm1;

    long int nextid = (myid + 1)%p;
    long int previd = (myid + npm1)%p;

    double complex psit1, psit2, psid1, psid2, psid3;
    for(long int lev = 1; lev<=log2np; lev++) {
        long int inc = p/pow(2, lev);
        for(long int j=1; j<=3; j++) {
            nextd[get1Dfrom2D(lev-1, j-1, p, 3)] = (myid + j*inc)%p;
            prevd[get1Dfrom2D(lev-1, j-1, p, 3)] = (myid + 3*p - j*inc)%p;
        }
    }

    for(long int box = 1; box<=t; box++) {
        long int ind = 1;
        long int lr = 1;

        for(long int sc=1; sc<=p; sc++) {
            for(long int j=1; j<=b; j++) {
                qrn[get1Dfrom2D(lr-1, ind-1, p, b)] = qr[get1Dfrom3D(box-1, sc-1, j-1, t, p, b)];
                lr += 1;
                if(lr>p) {
                    lr = 1;
                    ind = ind+1;
                }
            }
        }

        for(long int sc = 1; sc<=p; sc++) {
            for(long int j=1; j<=b; j++) {
                qr[get1Dfrom3D(box-1, sc-1, j-1, t, p, b)] = qrn[get1Dfrom2D(sc-1, j-1, p, b)];
            }
        }
    }

    // Get sum of each section.  Collect results in sumsec(1:npm1).
    for(long int sc=1; sc<=npm1; sc++) {
        sump[sc-1] = 0;

        for(long int box=1; box<=t; box++) {
            for(long int j=1; j<=b; j++) {
                sump[sc-1] = sump[sc-1] + qr[get1Dfrom3D(box-1, sc + 1 - 1, j-1, t, p, b)];
            }
        }
    }

    MPI_Reduce(sump, sumsec, npm1, MPI_C_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(sumsec, npm1, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    for(long int sc=1; sc<=npm1; sc++) {
        sumsec[sc-1] *= (I/dlq);
    }

    long int intrvl;
    long int intvlo;

    /* Step 1 : Lowest level */

    long int base = t + log2np - 3;
    initmm(qr, p, b, t, initch, &phi[get1Dfrom3D(base, 0, 0, 2*t + log2np-3, p-1, terms)], terms);


    /* Step 2 : Up the tree */

    long int lup = log2np;
    if(p==2)
        lup = 2;

    long int nl = t;
    //At lower levels, no communication required
    for(long int lev=n-1; lev>=lup; lev--) {
        long int obase = base;
        nl = nl/2;
        base = base -nl;

        for(long int box=1; box<=nl; box++) {
            mpp(&phi[get1Dfrom3D(obase + 2*box - 1 - 1, 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &shftfn[get1Dfrom3D(lev-1-1, 0, 0, n-2, terms, terms)], phils);
            mpp(&phi[get1Dfrom3D(obase + 2*box - 1, 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &shftfp[get1Dfrom3D(lev-1-1, 0, 0, n-2, terms, terms)], phirs);
            
            for(long int sc = 1; sc<=npm1;  sc++) {
                for(long int term = 1; term<=terms; term++) {
                    phi[get1Dfrom3D(base + box - 1, sc-1, term-1, 2*t + log2np - 3, p-1, terms)] = phils[get1Dfrom2D(sc-1, term-1, p-1, terms)] + phirs[get1Dfrom2D(sc-1, term-1, p-1, terms)];
                }
            }
        }
    }

    //At higher levels, communication is required
    for(long int lev = log2np - 1; lev>=2; lev--) {
        long int obase = base;
        base = base-1;
        intrvl = p/pow(2, lev);
        intvlo = intrvl/2;

        if(myid%intrvl == 0) {
            MPI_Recv(phr, pnpm1, MPI_C_DOUBLE_COMPLEX, nextd[get1Dfrom2D(lev+1-1, 0, p, 3)], lev, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            mpp(&phi[get1Dfrom3D(obase + 1 - 1 , 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &shftfn[get1Dfrom3D(lev-1-1, 0, 0, n-2, terms, terms)], phils);
            mpp(phr, terms, npm1, &shftfp[get1Dfrom3D(lev-1-1, 0, 0, n-2, terms, terms)], phirs);

            for(long int sc=1; sc<=npm1; sc++) {
                for(long int term=1; term<=terms; term++) {
                   phi[get1Dfrom3D(base + 1 - 1, sc-1, term-1, 2*t + log2np - 3, p-1, terms)] = phils[get1Dfrom2D(sc-1, term-1, p-1, terms)] + phirs[get1Dfrom2D(sc-1, term-1, p-1, terms)]; 
                }
            }
        }

        else if(myid % intrvl == intvlo) {
            MPI_Send(&phi[get1Dfrom3D(obase + 1 - 1 , 0, 0, 2*t + log2np - 3, p-1, terms)], pnpm1, MPI_C_DOUBLE_COMPLEX, prevd[get1Dfrom2D(lev+1-1, 0, p, 3)], lev, MPI_COMM_WORLD);
        }
    }

    long int tag = log2np;

    /* Pack:  phi, qr -> packns, packps */
    long int pk = 0;

    for(long int lev = lup+1; lev<=n; lev++) {
        long int nl = pow(2, lev)/p;
        base = nl + log2np - 3;

        for(long int sc=1; sc<=npm1; sc++) {
            for(long int term=1; term<=terms; term++) {
                pk++;
                packns[pk-1] = phi[get1Dfrom3D(base + nl - 1 - 1, sc-1, term-1, 2*t + log2np - 3, p-1, terms)];
                packps[pk-1] = phi[get1Dfrom3D(base + 1 - 1 , sc-1, term-1, 2*t + log2np - 3, p-1, terms)];

                pk++;
                packns[pk-1] = phi[get1Dfrom3D(base + nl - 1 , sc-1, term-1, 2*t + log2np - 3, p-1, terms)];
                packps[pk-1] = phi[get1Dfrom3D(base + 2 - 1 , sc-1, term-1, 2*t + log2np - 3, p-1, terms)];
            }
        }
    }

    /* Now pk .eq. (n - lup) * npm1 * 2 * p */

    for(long int sc=1; sc<=npm1; sc++) {
        for(long int j=1; j<=b; j++) {
            pk++;
            packns[pk-1] = qr[get1Dfrom3D(t-1, sc+1-1, j-1, t, p, b)];
            packps[pk-1] = qr[get1Dfrom3D(0, sc+1-1, j-1, t, p, b)];
        }
    }

    /* Now pk .eq. (n - lup)*npm1*2*p + npm1*nieach */
    long int lpkp = pk;

    if(t==1) {
        for(long int sc = p/2; sc<=npm1; sc++) {
            long int sce = sc - p/2 + 1;
            packnn[sce-1] = dotp(b, &qr[get1Dfrom3D(0, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 0, p/2, 3*b)]);
        }
    }
    else {
        for(long int sc=p/2; sc<=npm1; sc++) {
            long int sce = sc - p/2 + 1;
            pk++;
            packns[pk-1] = dotp(b, &qr[get1Dfrom3D(t-1-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 0, p/2, 3*b)]);
        }
    }

    long int lpkn = pk;

    /* Data transfer */
    tag++;
    MPI_Sendrecv(packps, lpkp, MPI_C_DOUBLE_COMPLEX, previd, tag, packnr, lpkp, MPI_C_DOUBLE_COMPLEX, nextid, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    tag++;
    MPI_Sendrecv(packns, lpkn, MPI_C_DOUBLE_COMPLEX, nextid, tag, packpr, lpkn, MPI_C_DOUBLE_COMPLEX, previd, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if(t==1) {
        tag++;
        MPI_Sendrecv(packnn, p/2, MPI_C_DOUBLE_COMPLEX, nextd[get1Dfrom2D(log2np-1, 1, p, 3)], tag, qcpp, p/2, MPI_C_DOUBLE_COMPLEX, prevd[get1Dfrom2D(log2np-1, 1, p, 3)], tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    /* Unpack:  packnr -> phin, qrn
          packpr -> phip, qrp, qcpp */

    pk = 0;
    long int lr = 0;

    for(long int lev = lup+1; lev<=n; lev++) {
        lr++;
        nl = pow(2, lev)/p;
        base = nl + log2np - 3;

        for(long int sc=1; sc<=npm1; sc++) {
            for(long int term=1; term<=terms; term++) {
                pk++;
                phin[get1Dfrom4D( lr-1, 0, sc-1, term-1, n, 2, p-1, terms)] = packnr[pk-1];
                phip[get1Dfrom4D( lr-1, 0, sc-1, term-1, n, 2, p-1, terms)] = packpr[pk-1];
                
                pk++;
                phin[get1Dfrom4D( lr-1, 1, sc-1, term-1, n, 2, p-1, terms)] = packnr[pk-1];
                phip[get1Dfrom4D( lr-1, 1, sc-1, term-1, n, 2, p-1, terms)] = packpr[pk-1];
            }
        }
    }

    /* Now pk .eq. (n - lup) * npm1 * 2 * p */
    for(long int sc=1; sc<=npm1; sc++) {
        for(long int j=1; j<=b; j++) {
            pk++;
            qrn[get1Dfrom2D(sc-1, j-1, p, b)] = packnr[pk-1];
            qrp[get1Dfrom2D(sc-1, j-1, p-1, b)] = packpr[pk-1];
        }
    }

    /* Now pk .eq. (n - lup)*npm1*2*p + nieach*npm1 */
    if(t>1) {
        for(long int sc=1; sc<=p/2; sc++) {
            pk++;
            qcpp[sc-1] = packpr[pk-1];
        }
    }
    
    /* Step 3-4:  Down the tree */
    /* -------- Top flip, which always requires communication */
    if(p==2) {
        tag++;
        nl = 2;
        for(long int box=1; box<=2; box++) {
            tag++;
            MPI_Sendrecv(&phi[get1Dfrom3D(box-1, 0, 0, 2*t + log2np - 3, p-1, terms)], pnpm1, MPI_C_DOUBLE_COMPLEX, nextid, tag, phiopp, pnpm1, MPI_C_DOUBLE_COMPLEX, previd, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            mpp(phiopp, terms, npm1, topflp, &psi[get1Dfrom3D(box-1, 0, 0, 2*t + log2np - 3, p-1, terms)]);
        }
    }
    else {
        tag++;
        nl = 1;
        intrvl = p/4;
        if(myid % intrvl == 0) {
            MPI_Sendrecv(&phi[get1Dfrom3D(0, 0, 0, 2*t + log2np - 3, p-1, terms)], pnpm1, MPI_C_DOUBLE_COMPLEX, nextd[get1Dfrom2D(0, 0, p, 3)], tag, phiopp, pnpm1, MPI_C_DOUBLE_COMPLEX, prevd[get1Dfrom2D(0, 0, p, 3)], tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            mpp(phiopp, terms, npm1, topflp, &psi[get1Dfrom3D(0, 0, 0, 2*t + log2np - 3, p-1, terms)]);
        }
    }

    base = 0;

 /* -------- Higher levels requiring communication */   
    for(long int lev = 3; lev<=log2np; lev++) {
        long int obase = base;
        base = base+1;
        intvlo = intrvl;
        intrvl = intvlo/2;
        tag+=10;

        if(myid%intvlo == 0) {
            MPI_Send(&psi[get1Dfrom3D(obase+1-1, 0, 0, 2*t + log2np - 3, p-1, terms)], pnpm1, MPI_C_DOUBLE_COMPLEX, nextd[get1Dfrom2D(lev-1, 0, p, 3)], tag+1, MPI_COMM_WORLD);
            MPI_Sendrecv(&phi[get1Dfrom3D(base+1-1, 0, 0, 2*t + log2np - 3, p-1, terms)], pnpm1, MPI_C_DOUBLE_COMPLEX, nextd[get1Dfrom2D(lev-1, 1, p, 3)], tag+2, phim2, pnpm1, MPI_C_DOUBLE_COMPLEX, prevd[get1Dfrom2D(lev-1, 1, p, 3)], tag+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&phi[get1Dfrom3D(base+1-1, 0, 0, 2*t + log2np - 3, p-1, terms)], pnpm1, MPI_C_DOUBLE_COMPLEX, prevd[get1Dfrom2D(lev-1, 1, p, 3)], tag+3, phip2, pnpm1, MPI_C_DOUBLE_COMPLEX, nextd[get1Dfrom2D(lev-1, 1, p, 3)], tag+3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&phi[get1Dfrom3D(base+1-1, 0, 0, 2*t + log2np - 3, p-1, terms)], pnpm1, MPI_C_DOUBLE_COMPLEX, nextd[get1Dfrom2D(lev-1, 2, p, 3)], tag+5, phip3, pnpm1, MPI_C_DOUBLE_COMPLEX, nextd[get1Dfrom2D(lev-1, 2, p, 3)], tag+4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            mpp(&psi[get1Dfrom3D(obase+1-1, 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &shftlp[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], shpsi);
            mpp(phim2, terms, npm1, &flip2p[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f2p);
            mpp(phip2, terms, npm1, &flip2n[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f2n);
            mpp(phip3, terms, npm1, &flip3n[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f3);

            for(long int sc=1; sc<=npm1; sc++) {
                for(long int term = 1; term<=terms; term++) {
                    psi[get1Dfrom3D(base+1-1, sc-1, term-1, 2*t + log2np - 3, p-1, terms)] = shpsi[get1Dfrom2D(sc-1, term-1, p-1, terms)] + f2p[get1Dfrom2D(sc-1, term-1, p-1, terms)] + f2n[get1Dfrom2D(sc-1, term-1, p-1, terms)] + f3[get1Dfrom2D(sc-1, term-1, p-1, terms)];
                }
            }
        }
        else if(myid%intvlo == intrvl) {
            MPI_Sendrecv(&phi[get1Dfrom3D(base+1-1, 0, 0, 2*t + log2np - 3, p-1, terms)], pnpm1, MPI_C_DOUBLE_COMPLEX, prevd[get1Dfrom2D(lev-1, 2, p, 3)], tag+4, psiprv, pnpm1, MPI_C_DOUBLE_COMPLEX, prevd[get1Dfrom2D(lev-1, 0, p, 3)], tag+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(phim3, pnpm1, MPI_C_DOUBLE_COMPLEX, prevd[get1Dfrom2D(lev-1, 2, p, 3)], tag+5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&phi[get1Dfrom3D(base+1-1, 0, 0, 2*t + log2np - 3, p-1, terms)], pnpm1, MPI_C_DOUBLE_COMPLEX, nextd[get1Dfrom2D(lev-1, 1, p, 3)], tag+6, phim2, pnpm1, MPI_C_DOUBLE_COMPLEX, prevd[get1Dfrom2D(lev-1, 1, p, 3)], tag+6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&phi[get1Dfrom3D(base+1-1, 0, 0, 2*t + log2np - 3, p-1, terms)], pnpm1, MPI_C_DOUBLE_COMPLEX, prevd[get1Dfrom2D(lev-1, 1, p, 3)], tag+7, phip2, pnpm1, MPI_C_DOUBLE_COMPLEX, nextd[get1Dfrom2D(lev-1, 1, p, 3)], tag+7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            mpp(psiprv, terms, npm1, &shftln[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], shpsi);
            mpp(phim2, terms, npm1, &flip2p[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f2p);
            mpp(phip2, terms, npm1, &flip2n[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f2n);
            mpp(phim3, terms, npm1, &flip3p[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f3);

            for(long int sc=1; sc<=npm1; sc++) {
                for(long int term=1; term<=terms; term++) {
                    psi[get1Dfrom3D(base+1-1, sc-1, term-1, 2*t + log2np - 3, p-1, terms)] = shpsi[get1Dfrom2D(sc-1, term-1, p-1, terms)] + f2p[get1Dfrom2D(sc-1, term-1, p-1, terms)] + f2n[get1Dfrom2D(sc-1, term-1, p-1, terms)] + f3[get1Dfrom2D(sc-1, term-1, p-1, terms)];
                }
            }
        }
    }

/* -------- Lower levels not requiring communication */

    lr = 0;
    for(long int lev=lup+1; lev<=n; lev++) {
        long int oldnl = nl;
        nl = nl*2;
        long int obase = base;
        base = base + oldnl;
        lr++;
        for(long int box=1; box<=oldnl; box++) {
            /* Left child */
            long int cl = base + 2*box - 1;  
            mpp(&psi[get1Dfrom3D(obase+box-1, 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &shftlp[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], shpsi);

            if(box==1) {
                mpp(&phip[get1Dfrom4D( lr-1, 0, 0, 0, n, 2, p-1, terms)], terms, npm1, &flip2p[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f2p);
            }
            else {
                mpp(&phi[get1Dfrom3D(cl-3, 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &flip2p[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f2p);
            }

            if(box==oldnl) {
                mpp(&phin[get1Dfrom4D( lr-1, 0, 0, 0, n, 2, p-1, terms)], terms, npm1, &flip2n[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f2n);
                mpp(&phin[get1Dfrom4D( lr-1, 1, 0, 0, n, 2, p-1, terms)], terms, npm1, &flip3n[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f3);
            }
            else {
                mpp(&phi[get1Dfrom3D(cl+1, 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &flip2n[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f2n);
                mpp(&phi[get1Dfrom3D(cl+2, 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &flip3n[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f3);
            }

            for(long int sc=1; sc<=npm1; sc++) {
                for(long int term=1; term<=terms; term++) {
                    psi[get1Dfrom3D(cl-1, sc-1, term-1, 2*t + log2np - 3, p-1, terms)] = shpsi[get1Dfrom2D(sc-1, term-1, p-1, terms)] + f2p[get1Dfrom2D(sc-1, term-1, p-1, terms)] + f2n[get1Dfrom2D(sc-1, term-1, p-1, terms)] + f3[get1Dfrom2D(sc-1, term-1, p-1, terms)];
                }
            }

            long int cr = base + 2* box;
            
            mpp(&psi[get1Dfrom3D(obase+box-1, 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &shftln[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], shpsi);

            if(box==1) {
                mpp(&phip[get1Dfrom4D( lr-1, 0, 0, 0, n, 2, p-1, terms)], terms, npm1, &flip3p[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f3);
                mpp(&phip[get1Dfrom4D( lr-1, 1, 0, 0, n, 2, p-1, terms)], terms, npm1, &flip2p[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f2p);
            }
            else {
                mpp(&phi[get1Dfrom3D(cr-4, 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &flip3p[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f3);
                mpp(&phi[get1Dfrom3D(cr-3, 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &flip2p[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f2p);
            }

            if(box==oldnl) {
                mpp(&phin[get1Dfrom4D( lr-1, 1, 0, 0, n, 2, p-1, terms)], terms, npm1, &flip2n[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f2n);
            }
            else {
                mpp(&phi[get1Dfrom3D(cr+1, 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &flip2n[get1Dfrom3D(lev-3, 0, 0, n-2, terms, terms)], f2n);
            }

            for(long int sc=1; sc<=npm1; sc++) {
                for(long int term=1; term<=terms; term++) {
                    psi[get1Dfrom3D(cr-1, sc-1, term-1, 2*t + log2np - 3, p-1, terms)] = shpsi[get1Dfrom2D(sc-1, term-1, p-1, terms)] + f2p[get1Dfrom2D(sc-1, term-1, p-1, terms)] + f2n[get1Dfrom2D(sc-1, term-1, p-1, terms)] + f3[get1Dfrom2D(sc-1, term-1, p-1, terms)];
                }
            }
        }
    }

    /* Step 5:  Evaluate local expansions. */
    
    for(long int sc=p/2; sc<=npm1; sc++) {
        long int sce = sc - p/2 + 1;
        exts[sce-1] = dotp(terms, &psi[get1Dfrom3D(base + t - 1, sc-1, 0, 2*t + log2np - 3, p-1, terms)], &evalmh[get1Dfrom2D(sce-1, 0, p/2, terms)]);   
    }

    tag+=10;

    MPI_Send(exts, p/2, MPI_C_DOUBLE_COMPLEX, nextid, tag, MPI_COMM_WORLD);
    MPI_Recv(extr, p/2, MPI_C_DOUBLE_COMPLEX, previd, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for(long int sc=1; sc<=p/2 -1; sc++) {
        for(long int j=1; j<=b; j++) {
            for(long int box=1; box<=t; box++) {
                psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] = dotp(terms, &psi[get1Dfrom3D(base+box-1, sc-1, 0, 2*t + log2np - 3, p-1, terms)], &evalm[get1Dfrom3D(j-1, sc-1, 0, b, p-1, terms)]);
            }
        }
    }

    long int sc = p/2;
    long int sce = 1;
    long int j=1;
    long int box = 1;
    psit1 = dotp(terms, &psi[get1Dfrom3D(base+box-1, sc-1, 0, 2*t + log2np - 3, p-1, terms)], &evalm[get1Dfrom3D(j-1, sc-1, 0, b, p-1, terms)]);
    psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] = 0.5 * (psit1 + extr[sce-1]);

    for(box = 2; box<=t; box++) {
        psit1 = dotp(terms, &psi[get1Dfrom3D(base+box-1, sc-1, 0, 2*t + log2np - 3, p-1, terms)], &evalm[get1Dfrom3D(j-1, sc-1, 0, b, p-1, terms)]);
        psit2 = dotp(terms, &psi[get1Dfrom3D(base+box-2, sc-1, 0, 2*t + log2np - 3, p-1, terms)], &evalmh[get1Dfrom2D(sce-1, 0, p/2, terms)]);
        psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] = 0.5 * (psit1 + psit2);
    }

    for(j=2; j<=b; j++) {
        for(box=1; box<=t; box++) {
            psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] = dotp(terms, &psi[get1Dfrom3D(base+box-1, sc-1, 0, 2*t + log2np - 3, p-1, terms)], &evalm[get1Dfrom3D(j-1, sc-1, 0, b, p-1, terms)]);
        }
    }

    for(sc=p/2+1; sc<=npm1; sc++) {
        sce = sc - p/2 + 1;
        j=1;
        box = 1;
        psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] = extr[sce-1];
        for(box=2; box<=t; box++) {
            psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] = dotp(terms, &psi[get1Dfrom3D(base+box-2, sc-1, 0, 2*t + log2np - 3, p-1, terms)], &evalmh[get1Dfrom2D(sce-1, 0, p/2, terms)]);
        }
        for(j=2; j<=b; j++) {
            for(box = 1; box<=t; box++) {
                psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] = dotp(terms, &psi[get1Dfrom3D(base+box-1, sc-1, 0, 2*t + log2np - 3, p-1, terms)], &evalm[get1Dfrom3D(j-1, sc-1, 0, b, p-1, terms)]);
            }
        }
    }

    /* Step 6:  Direct evaluation. */
    for(long int sc=1; sc<=p/2-1; sc++) {
        j=1;
        if(t==1) {
            box=1;
            psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qrn[get1Dfrom2D(sc-1,0,p, b)], &cotprv[get1Dfrom2D(sc-1, 0, p/2, 3*b)]);
        }
        else {
            box = 1;
            psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box+1-1, sc+1-1, 0, t, p, b)], &cotprv[get1Dfrom2D(sc-1, 0, p/2, 3*b)]);
            for(box=2; box<=t-1; box++) {
               psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box+1-1, sc+1-1, 0, t, p, b)], &cotprv[get1Dfrom2D(sc-1, 0, p/2, 3*b)]);
            }
            box = t;
            psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qrn[get1Dfrom2D(sc-1,0,p, b)], &cotprv[get1Dfrom2D(sc-1, 0, p/2, 3*b)]);
        }

        for(long int j=2; j<=b; j++) {
            if(t==1) {
                box=1;
                psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qrn[get1Dfrom2D(sc-1,0,p, b)], &cots[get1Dfrom3D(sc-1,j-2,0,p-1, b-1, 3*b)]);
            }
            else {
                box=1;
                psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box+1-1, sc+1-1, 0, t, p, b)], &cots[get1Dfrom3D(sc-1,j-2,0,p-1, b-1, 3*b)]);
                for(box=2; box<=t-1; box++) {
                   psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box+1-1, sc+1-1, 0, t, p, b)], &cots[get1Dfrom3D(sc-1,j-2,0,p-1, b-1, 3*b)]); 
                }
                box = t;
                psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qrn[get1Dfrom2D(sc-1,0,p, b)], &cots[get1Dfrom3D(sc-1,j-2,0,p-1, b-1, 3*b)]);
            }
        }
    }

    sc = p/2;
    sce = 1;
    j = 1;
    if(t==1) {
        box=1;
        psid1 = mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qrn[get1Dfrom2D(sc-1,0,p, b)], &cotprv[get1Dfrom2D(sc-1, 0, p/2, 3*b)]);
        psid2 = dotp(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &cotsh[get1Dfrom2D(sce-1, 1+b-1, p/2, 3*b)]);
        psid3 = dotp(b, &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 1+2*b-1, p/2, 3*b)]);
        psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += (psid1 + qcpp[sce-1] + psid2 + psid3) * 0.5;
    }
    else if(t==2) {
        box=1;
        psid1 = mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box+1-1, sc+1-1, 0, t, p, b)], &cotprv[get1Dfrom2D(sc-1, 0, p/2, 3*b)]);
        psid2 = dotp(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &cotsh[get1Dfrom2D(sce-1, 1+b-1, p/2, 3*b)]);
        psid3 = dotp(b,&qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 1+2*b-1, p/2, 3*b)]);
        psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += (psid1 + qcpp[sce-1] + psid2 + psid3) * 0.5;

        box=2;
        psid1 = mn3(b, &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qrn[get1Dfrom2D(sc-1,0,p, b)], &cotprv[get1Dfrom2D(sc-1, 0, p/2, 3*b)]);
        psid2 = mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 0, p/2, 3*b)]);
        psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += (psid1 + psid2) * 0.5;
    }
    else {
        box = 1;
        psid1 = mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box+1-1, sc+1-1, 0, t, p, b)], &cotprv[get1Dfrom2D(sc-1, 0, p/2, 3*b)]);
        psid2 = dotp(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &cotsh[get1Dfrom2D(sce-1, 1+b-1, p/2, 3*b)]);
        psid3 = dotp(b,&qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 1+2*b-1, p/2, 3*b)]);
        psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += (psid1 + qcpp[sce-1] + psid2 + psid3) * 0.5;

        box = 2;
        psid1 = mn3(b, &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box+1-1, sc+1-1, 0, t, p, b)], &cotprv[get1Dfrom2D(sc-1, 0, p/2, 3*b)]);
        psid2 = mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 0, p/2, 3*b)]);
        psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += (psid1 + psid2) * 0.5;

        for(long int box=3; box<=t-1; box++) {
          psid1 = mn3(b, &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box+1-1, sc+1-1, 0, t, p, b)], &cotprv[get1Dfrom2D(sc-1, 0, p/2, 3*b)]);
          psid2 = mn3(b, &qr[get1Dfrom3D(box-3, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 0, p/2, 3*b)]);
          psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += (psid1 + psid2) * 0.5;
        }

        box = t;
        psid1 = mn3(b, &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qrn[get1Dfrom2D(sc-1,0,p, b)], &cotprv[get1Dfrom2D(sc-1, 0, p/2, 3*b)]);
        psid2 = mn3(b, &qr[get1Dfrom3D(box-3, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 0, p/2, 3*b)]);
        psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += (psid1 + psid2) * 0.5;
    }

    for(j=2; j<=b; j++) {
        if(t==1) {
            box = 1;
            psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qrn[get1Dfrom2D(sc-1,0,p, b)], &cots[get1Dfrom3D(sc-1,j-2,0,p-1, b-1, 3*b)]);
        }
        else {
            box = 1;
            psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box, sc+1-1, 0, t, p, b)], &cots[get1Dfrom3D(sc-1,j-2,0,p-1, b-1, 3*b)]);
            for(box = 2; box<=t-1; box++) {
                psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box, sc+1-1, 0, t, p, b)], &cots[get1Dfrom3D(sc-1,j-2,0,p-1, b-1, 3*b)]);
            }
            box = t;
            psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qrn[get1Dfrom2D(sc-1,0,p, b)], &cots[get1Dfrom3D(sc-1,j-2,0,p-1, b-1, 3*b)]);
        }
    }

    for(sc = p/2+1; sc<=npm1; sc++) {
        sce = sc - p/2 + 1;
        j=1;

        if(t==1) {
            box = 1;
            psid1 = dotp(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &cotsh[get1Dfrom2D(sce-1, b, p/2, 3*b)]);
            psid2 = dotp(b, &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 2*b, p/2, 3*b)]);
            psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += qcpp[sce-1] + psid1 + psid2;
        }
        else if(t==2) {
            box = 1;
            psid1 = dotp(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &cotsh[get1Dfrom2D(sce-1, b, p/2, 3*b)]);
            psid2 = dotp(b, &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 2*b, p/2, 3*b)]);
            psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += qcpp[sce-1] + psid1 + psid2;

            box = 2;
            psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 0, p/2, 3*b)]);
        }
        else {
            box = 1;
            psid1 = dotp(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &cotsh[get1Dfrom2D(sce-1, b, p/2, 3*b)]);
            psid2 = dotp(b, &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 2*b, p/2, 3*b)]);
            psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += qcpp[sce-1] + psid1 + psid2;

            box = 2;
            psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 0, p/2, 3*b)]);

            for(long int box = 3; box<=t-1; box++) {
                psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qr[get1Dfrom3D(box-3, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 0, p/2, 3*b)]);
            }

            box = t;
            psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qr[get1Dfrom3D(box-3, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 0, p/2, 3*b)]);
        }

        for(j=2; j<=b; j++) {
            if(t==1) {
                box=1;
                psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qrn[get1Dfrom2D(sc-1,0,p, b)], &cots[get1Dfrom3D(sc-1,j-2,0,p-1, b-1, 3*b)]);
            }
            else {
                box = 1;
                psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qrp[get1Dfrom2D(sc-1, 0, p-1, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box, sc+1-1, 0, t, p, b)], &cots[get1Dfrom3D(sc-1,j-2,0,p-1, b-1, 3*b)]);

                for(long int box = 2; box<=t-1; box++) {
                  psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box, sc+1-1, 0, t, p, b)], &cots[get1Dfrom3D(sc-1,j-2,0,p-1, b-1, 3*b)]);
                }

                box = t;
                psiev[get1Dfrom3D(sc-1, box-1, j-1, p, t, b)] += mn3(b, &qr[get1Dfrom3D(box-2, sc+1-1, 0, t, p, b)], &qr[get1Dfrom3D(box-1, sc+1-1, 0, t, p, b)], &qrn[get1Dfrom2D(sc-1,0,p, b)], &cots[get1Dfrom3D(sc-1,j-2,0,p-1, b-1, 3*b)]);
            }
        }
    }

    for(long int box=1; box<=t; box++) {
        for(long int pk=1; pk<=b; pk++) {
            for(long int j=1; j<=npm1; j++) {
                psit[j-1] = psiev[get1Dfrom3D(j-1, box-1, pk-1, p, t, b)];
            }
            for(long int j=1; j<=p; j++) {
                psiev[get1Dfrom3D(j-1, box-1, pk-1, p, t, b)] = qr[get1Dfrom3D(box-1, 0, pk-1, t, p, b)];
                for(long int sc=1; sc<=npm1; sc++) {
                    psiev[get1Dfrom3D(j-1, box-1, pk-1, p, t, b)] += fo[get1Dfrom2D(j-1, sc-1, p, p-1)] * (psit[sc-1] - sumsec[sc-1]);
                }
            }
        }
    }

    MPI_Alltoall(psiev, t*b, MPI_C_DOUBLE_COMPLEX, qr, t*b, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    //myfft();

    if(dir<0)
        cjall(lq, qr);
}

void mft(long int lq, double complex* qr, long int dir, long int p, long int myid, long int terms, long int b, double* w, double* v, long int sz1, long int szwk, fftw_plan forward_plan) {
    MPI_Barrier(MPI_COMM_WORLD);
    long int s = lq/b;
    long int t = s/p;
    long int n = my_log2(s);
    long int log2np = my_log2(p);

    long int ia[57];
    long int ind = 0;

    for(long int j=0; j<=7; j++) {
        ia[j] = ind;
        ind += terms*terms*(n-2);
    }
    ia[8] = ind;
    ind += terms*terms;
    ia[9] = ind;
    ind += terms*b;
    ia[10] = ind;
    ind += terms*(p - 1)*b;
    ia[11] = ind;
    ind += terms*p/2;
    ia[12] = ind;
    ind += (3*b*p)/2;
    ia[13] = ind;
    ind += 3*b*(b-1)*(p-1);
    ia[14] = ind;
    ind += (3*b*p)/2;
    ia[15] = ind;
    ind += terms;
    ia[16] = ind;
    ind += 2*(p-1)*p;
    ind = 0;
    for(long int j=17; j<=18; j++) {
        ia[j] = ind;
        ind+= 2*terms*(p-1)*(2*t + log2np - 3);
    }
    for(long int j=21; j<=22; j++) {
        ia[j] = ind;
        ind += 2*p;
    }
    for(long int j=23; j<=25; j++) {
        ia[j] = ind;
        ind += 2*terms*(p-1);
    }
    for(long int j=26; j<=27; j++) {
        ia[j] = ind;
        ind += 2*(p-1)*(b + n*2*terms);
    }
    for(long int j=28; j<=29; j++) {
        ia[j] = ind;
        ind += 2*((p-1)*(b + n*2*terms) + p/2);
    }
    ia[30] = ind;
    ind += 2*b*p;
    ia[31] = ind;
    ind += 2*b*(p-1);
    for(long int j=32; j<=33; j++) {
        ia[j] = ind;
        ind+= 2 * (p/2);
    }
    for(long int j=34; j<=43; j++) {
        ia[j] = ind;
        ind += 2*terms*(p-1);
    }
    for(long int j=44; j<=45; j++) {
        ia[j] = ind;
        ind += 2 * terms * (p-1) * 2 * n;
    }
    ia[46] = ind;
    ind += 2*lq;
    ia[47] = ind;
    ind += 2*(p-1);
    for(long int j=48; j<=49; j++) {
        ia[j] = ind;
        ind += 2 * (p/2);
    }
    ia[50] = ind;
    ind = ind + b;
    for(long int j=51; j<=53; j++) {
        ia[j] = ind;
        ind += terms;
    }
    for(long int j=54; j<=55; j++) {
        ia[j] = ind;
        ind += 3 * p;
    }

    mftint(qr, lq, p, myid, s, terms, n, t, b, log2np, sz1, szwk, dir,
        &w[ia[0]], &w[ia[1]], &w[ia[2]], &w[ia[3]], &w[ia[4]], &w[ia[5]],
        &w[ia[6]], &w[ia[7]], &w[ia[8]], &w[ia[9]], &w[ia[10]], &w[ia[11]], 
        &w[ia[12]], &w[ia[13]], &w[ia[14]], &w[ia[15]], &w[ia[16]],
        &v[ia[17]], &v[ia[18]], &v[ia[21]], &v[ia[22]], &v[ia[23]], &v[ia[24]],
        &v[ia[25]], &v[ia[26]], &v[ia[27]], &v[ia[28]], &v[ia[29]], &v[ia[30]],
        &v[ia[31]], &v[ia[32]], &v[ia[33]], &v[ia[34]], &v[ia[35]], &v[ia[36]],
        &v[ia[37]], &v[ia[38]], &v[ia[39]], &v[ia[40]], &v[ia[41]], &v[ia[42]],
        &v[ia[43]], &v[ia[44]], &v[ia[45]], &v[ia[46]], &v[ia[47]], &v[ia[48]],
        &v[ia[49]], &v[ia[50]], &v[ia[51]], &v[ia[52]], &v[ia[53]], (void*) &v[ia[54]], (void*) &v[ia[55]]
    );

    fftw_execute(forward_plan);
}


void display(fftw_complex* x, fftw_complex* y, long int len) {
    for(long int i=0; i<len; i++) {
        /* term wise error */
        double err = cabs(x[i]-y[i])/cabs(y[i]) * 100;

        printf("[%d] (%.1f, %.1f) ---- (%.1f, %.1f) \t\t %lf\n", 
            i, creal(x[i]), cimag(x[i]), 
            creal(y[i]), cimag(y[i]), 
            err
        );
    }
}

double maxError(fftw_complex* x, fftw_complex* y, long int len) {
    double maxErr = cabs(x[0]-y[0])/cabs(y[0]) * 100;

    for(long int i=1; i<len; i++) {
        /* term wise error */
        double err = cabs(x[i]-y[i])/cabs(y[i]) * 100;
        if(maxErr < err) {
            maxErr = err;
        }
    }
    return maxErr;
}

long int getvlen(long int lq, long int terms, long int p, long int b) {
    long int s = lq/b;
    long int t = s/p;
    long int n = my_log2(s);
    long int log2np = my_log2(p);

    long int ind = 0;
    ind += 2*(2*terms*(p-1)*(2*t + log2np - 3));
    ind += 2*(2*p);
    ind += 3*(2*terms*(p-1));
    ind += 2*(2*(p-1)*(b + n*2*terms));
    ind += 2*(2*((p-1)*(b + n*2*terms) + p/2));
    ind += 2*b*p;
    ind += 2*b*(p-1);
    ind += 2* (2 * (p/2));
    ind += 10*(2*terms*(p-1));
    ind += 2*(2 * terms * (p-1) * 2 * n);
    ind += 2*lq;
    ind += 2*(p-1);
    ind += 2*(2 * (p/2));
    ind = ind + b;
    ind += 3*terms;
    ind += 2 * (3 * p);
    ind++;
    return ind;
}

long int getwlen(long int lq, long int terms, long int p, long int b){
    long int s = lq/b;
    long int t = s/p;
    long int n = my_log2(s);
    long int log2np = my_log2(p);

    long int ind = 0;
    ind += 8*(terms*terms*(n-2));
    ind += terms*terms;
    ind += terms*b;
    ind += terms*(p - 1)*b;
    ind += terms*p/2;
    ind += (3*b*p)/2;
    ind += 3*b*(b-1)*(p-1);
    ind += (3*b*p)/2;
    ind += terms;
    ind += 2*(p-1)*p;
    ind += p*(p-1);
    ind++;
    
    return ind;
}

long int getwlen2(long int lq, long int terms, long int p, long int b){
    long int s = lq/b;
    long int t = s/p;
    long int n = my_log2(s);
    long int log2np = my_log2(p);

    long int ind = 0;
    ind += 8*(terms*terms*(n-2));
    ind += terms*terms;
    ind += terms*b;
    ind += terms*(p - 1)*b;
    ind += terms*p/2;
    ind += (3*b*p)/2;
    ind += 3*b*(b-1)*(p-1);
    ind += (3*b*p)/2;
    ind += terms;
    ind += 2*(p-1)*p;
    ind += p*(p-1);
    ind+=b;
    ind+= 3*terms;
    ind++;
    
    return ind;
}


/*
* N : fft size
* P : number of processors
* B : number of boxes
* T : number of terms
*/
long int main(long int argc, char* argv[]) {
    long int myid, world_size;
    long int N, P, B, T;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);

    if(argc<3) {
        printf("ERROR : Provide correct args\n");
        MPI_Finalize();
        return 0;
    }

    N = atoi(argv[1]);
    B = atoi(argv[2]);
    T = atoi(argv[3]);

    P = world_size;

    long int local_length = N/P;

    fftw_complex *x = 0;
    x  = fftw_malloc(sizeof(fftw_complex)*local_length);

    fftw_plan forward_plan = fftw_plan_dft(1, &local_length, x, x, FFTW_FORWARD, FFTW_ESTIMATE);

    long int w_elements, v_elements;
    double* w;
    double* v;
    w_elements = getwlen2(N/P, T, P, B);
    v_elements = getvlen(N/P, T, P, B);
    w = (double*)fftw_malloc(sizeof(double) * w_elements);
    v = (double*)fftw_malloc(sizeof(double) * v_elements);

    /* initialise signal */
    if(myid == 0) {
        /* send */
        long int count = 0;
        for(long int proc = 1; proc<P; proc++) {
            for(long int i=0; i<local_length; i++) {
                x[i] = (fftw_complex)(proc * local_length + i + 1);
                count++;
            }
            MPI_Send(x, local_length, MPI_C_DOUBLE_COMPLEX, proc, 0, MPI_COMM_WORLD);
        }
        /* assign local data of id 0 */
        for(long int i=0; i<local_length; i++) {
            x[i] = i+1;
        }
    }
    else {
        /* receive */
        MPI_Recv(x, local_length, MPI_C_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    /* mfti */
    mfti(local_length, P, T, B, w, v, 0, 0);

    long int RUNS = 6;
    double global_time = 0;

    for(long int runid = 0; runid<RUNS; runid++) {
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* mft */
        // double start = MPI_Wtime();
        double start = omp_get_wtime();
        mft(local_length, x, 1, P, myid, T, B, w, v, 0, 0, forward_plan);
        MPI_Barrier(MPI_COMM_WORLD);
        // double end = MPI_Wtime() - start;
        double end = omp_get_wtime() - start; 

        if(myid==0) {
            printf("%d, %d, %d, %d, %lf, %lf\n", N, P, B, T, end, 0);
        }

        if(runid>0)
            global_time+=end;
        
        MPI_Barrier(MPI_COMM_WORLD); 
    }
    
    global_time = global_time/(RUNS-1);

    // if(myid==0) {
    //     printf("%d, %d, %d, %d, %lf, %lf\n", N, P, B, T, global_time, 0);
    // }

    /* free resources */
    MPI_Barrier(MPI_COMM_WORLD);
    fftw_destroy_plan(forward_plan);
    fftw_free(x);
    fftw_free(w);
    fftw_free(v);
    MPI_Finalize();
}
