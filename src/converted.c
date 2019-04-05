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

int get1Dfrom4D(int i, int j, int k, int l, int dim1, int dim2, int dim3, int dim4) {
    return((i * dim2 * dim3 * dim4) + (j*dim3*dim4) + k*dim4 + l);
}

int my_log(int x) {
    int ans = 0;
    while(x>1) {
        x/=2;
        ans++;
    }
    return ans;
}

void cjall(int n, double complex* x) {
    for(int i=0; i<n; i++) {
        x[i] = conj(x[i]);
    }
}

double complex dotp(int n, double complex* x, double* r) {
    double complex ans = 0;
    for(int i=0; i<n; i++) {
        ans += x[i]*r[i];
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

void wx(double a, double* x, int xlen, double* mat, int terms, double* chnods, double* acu) {
    double ac = 1/tan(a);

    for(int j=1; j<=terms; j++) {
        double th = -cos( ((double)(2*j - 2)) / 2*(terms-1) * M_PI ) * (M_PI - 6*a) + M_PI;
        acu[j-1] = ac * tan(0.5 * th);
        chnods[j-1] = -cos(((double)2*j -1) / (2*terms) * M_PI);
    }

    for(int j=1; j<=terms; j++) {
        double mm = chnods[j-1] - acu[j-1];

        for(int k=1; k<=j-1; k++) {
            mm *= ( (chnods[j-1] - acu[k-1]) / (chnods[j-1] - chnods[k-1]) );
        }

        for(int k=j+1; k<=terms; k++) {
            mm *= ( (chnods[j-1] - acu[k-1]) / (chnods[j-1] - chnods[k-1]) );
        }
        
        for(int k=1; k<=xlen; k++) {
            mat[get1Dfrom2D(k-1, j-1, xlen, terms)] = mm / (x[k-1] - acu[j-1]);
            for(int l=1; l<=j-1; l++) {
                mat[get1Dfrom2D(k-1, j-1, xlen, terms)] *= ( (x[k-1] - chnods[l-1]) / (x[k-1] - acu[l-1]) );
            }
            for(int l=j+1; l<=terms; l++) {
                mat[get1Dfrom2D(k-1, j-1, xlen, terms)] *= ( (x[k-1] - chnods[l-1]) / (x[k-1] - acu[l-1]) );
            }
        }
    }
}

void shftf(int lev, int dir, double* mat, int terms, double* chnods, double* x, double* wkp, double* mm) {
    double dd = (double) dir;
    double a = M_PI/pow(2, lev);
    double ta3 = 2 * tan(a);
    double td2 = tan(0.5 * a);

    for(int j=1; j<=terms; j++) {
        x[j-1] = 3 * td2 * (chnods[j-1] + dd*ta3*td2) / (ta3 - dd*td2*chnods[j-1]);
    }

    vx(0.5*a, x, mat, terms, chnods, wkp, mm);
}

void evlmf(double ffr, int s, int lq, double* vec, int terms, double* w, double* chnods, double* acu) {
    double a = (M_PI / s) * 0.5;
    double th = ffr * a;
    double x = tan(th)/tan(a);
    wx(a, &x, 1, w, terms, chnods, acu);

    for(int t = 1; t<=terms; t++) {
        vec[t-1] = w[t-1]/lq;
    }

    return;
}

void mpp(double complex* vecs, int terms, int len, double* matpp, double complex* ans) {
    for(int sc=1; sc<=len; sc++) {
        for(int term=1; term<=terms; term++) {
            ans[get1Dfrom2D(sc-1, term-1,len, terms)] = vecs[get1Dfrom2D(sc-1, 0,len ,terms)] * matpp[get1Dfrom2D(term-1, 0, terms, terms)];

            for(int k=2; k<=terms; k++) {
                ans[get1Dfrom2D(sc-1, term-1,len, terms)] += vecs[get1Dfrom2D(sc-1, k-1,len ,terms)] * matpp[get1Dfrom2D(term-1, k-1, terms, terms)];
            }
        }
    }
}

void shftl(int lev, int dir, double* mat, int terms, double* chnods, double* x, double* acu) {
    double dd = (double)dir;
    double a = M_PI/pow(2, lev);
    double c = tan(0.25 * a);
    double td2 = tan(0.5 * a);

    for(int j=1; j<=terms; j++) {
        x[j-1] = (chnods[j] - dd) / (1/c + dd*c*chnods[j-1]) / td2;
    }

    wx(0.5 * a, x, terms, mat, terms, chnods, acu);
}

void initmm(double complex* qr, int p, int b, int t, double* initch, double complex* phi, int terms) {
    for(int box=1; box<=t; box++) {
        for(int sc=1; sc<=p-1; sc++) {
            for(int term=1; term<=terms; term++) {
                phi[get1Dfrom3D(box-1, sc-1, term-1, t, p-1, terms)] = qr[get1Dfrom3D(box-1, sc, 0, t, p, b)] * initch[get1Dfrom2D(0, term-1, b, terms)];

                for(int j=2; j<=b; j++) {
                    phi[get1Dfrom3D(box-1, sc-1, term-1, t, p-1, terms)] += qr[get1Dfrom3D(box-1, sc, j-1, t, p, b)] * initch[get1Dfrom2D(j-1, term-1, b, terms)];
                }
            }
        }
    }
}

void flip(int lev, int shft, double* mat, int terms, double* chnods, double* x, double* wkp, double* mm) {
    double a = M_PI/(pow(2, lev));
    double b = shft;
    double c = tan(0.5 * a);
    double td2 = tan(b * a);

    for(int j=1; j<=terms; j++) {
        x[j-1] = 3 * c * (1 - td2*c*chnods[j-1]) / (td2 + c*chnods[j-1]);
    }

    vx(0.5 * a, x, mat, terms, chnods, mm, wkp);

    return;
}

void initg(int b, int s, double* initch, int terms, double* tang, double* chnods) {
    double a = M_PI / (2.0 * s);

    for(int j=1; j<=b; j++) {
        double th = ((double)(2*j - b - 1))/(double)b * a;
        tang[j-1] = tan(th);
    }

    double ta3 = 3 * tan(a);

    for(int j=1; j<=b; j++) {
        for(int k=1; k<=terms; k++) {
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

        for(int j=1; j<=terms; j++) {
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
                evlmf(ffr, s, lq, &evalm[get1Dfrom3D(j-1, sc-1, 0, b, p-1, terms)], terms, wkp, chnods, wkp2);
            }

            if(sc >= p/2) {
                int sce = sc - p/2 + 1;
                double ffr = (((double)2*b + 1) - ((double)2*sc)/(p))/(b) - 1;
                evlmf(ffr, s, lq, &evalmh[get1Dfrom2D(sce-1, 0, p/2, terms)], terms, wkp, chnods, wkp2);

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
                double complex temp = (ang * (1 - 2*j))*I;
                fo[get1Dfrom2D(sc-1, j+b-1, p, p-1)] = -cexp(temp)*sin(ang);
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

void mftint(double complex* qr, int lq, int p, int myid, int s, int terms, int n, int t, int b, int log2np, int sz1, int szwk, int dir,
    double* shftfn, double* shftfp, double* flip2n, double* flip2p, double* flip3n, double* flip3p,
    double* shftln, double* shftlp, double* topflp, double* initch, double* evalm, double* evalmh,
    double* cotsh, double* cots, double* cotprv, double* cnods, double complex* fo,
    double complex* phi, double complex* psi, double complex* sump, double complex* sumsec, double complex* phils, double complex* phirs, double complex* phr,
    double complex* packps, double complex* packnr, double complex* packns, double complex* packpr, double complex* qrn, double complex* qrp, double complex* qcpp, double complex* packnn,
    double complex* phiopp, double complex* phim2, double complex* phim3, double complex* phip3, double complex* shpsi, double complex* psiprv,
    double complex* f2n, double complex* f2p, double complex* f3, double complex* phin, double complex* phip, double complex* psiev, double complex* psit, double complex* exts, double complex* extr,
    double* wknie, double* wkp, double* wkp2, double* wkp3,
    int* prevd, int* nextd
) {

    MPI_Barrier(MPI_COMM_WORLD);

    if(dir<0) {
        cjall(lq, qr);
    }

    double dlq = (double)lq;
    log2np = my_log2(p);
    int npm1 = p-1;
    int pnpm1 = terms*npm1;

    int nextid = (myid + 1)%p;
    int previd = (myid - 1)%p;

    for(int lev = 1; lev<=log2np; lev++) {
        int inc = p/pow(2, lev);
        for(int j=1; j<=3; j++) {
            nextd[get1Dfrom2D(lev-1, j-1, p, 3)] = (myid + j*inc)%p;
            prevd[get1Dfrom2D(lev-1, j-1, p, 3)] = (myid + 3*p - j*inc)%p;
        }
    }

    for(int box = 1; box<=t; box++) {
        int ind = 1;
        int lr = 1;

        for(int sc=1; sc<=p; sc++) {
            for(int j=1; j<=b; j++) {
                qrn[get1Dfrom2D(lr-1, ind-1, p, b)] = qr[get1Dfrom3D(box-1, sc-1, j-1, t, p, b)];
                lr += 1;
                if(lr>p) {
                    lr = 1;
                    ind = ind+1;
                }
            }
        }

        for(int sc = 1; sc<=p; sc++) {
            for(int j=1; j<=b; j++) {
                qr[get1Dfrom3D(box-1, sc-1, j-1, t, p, b)] = qrn[get1Dfrom2D(sc-1, j-1, p, b)];
            }
        }
    }

    // Get sum of each section.  Collect results in sumsec(1:npm1).
    for(int sc=1; sc<=npm1; sc++) {
        sump[sc-1] = 0;

        for(int box=1; box<=t; box++) {
            for(int j=1; j<=b; j++) {
                sump[sc-1] = sump[sc-1] + qr[get1Dfrom3D(box-1, sc + 1 - 1, j-1, t, p, b)];
            }
        }
    }

    MPI_Reduce(sump, sumsec, npm1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(sumsec, npm1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    for(int sc=1; sc<=npm1; sc++) {
        sumsec[sc-1] *= (I/dlq);
    }

    /* Step 1 : Lowest level */

    int base = t + log2np - 3;
    initmm(qr, p, b, t, initch, phi, terms);

    /* Step 2 : Up the tree */

    int lup = log2np;
    if(p==2)
        lup = 2;

    int nl = t;
    //At lower levels, no communication required
    for(int lev=n-1; lev>=lup; lev--) {
        int obase = base;
        nl = nl/2;
        base = base -nl;

        for(int box=1; box<=nl; box++) {
            mpp(&phi[get1Dfrom3D(obase + 2*box - 1 - 1, 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &shftfn[get1Dfrom3D(lev-1-1, 0, 0, n-2, terms, terms)], phils);
            mpp(&phi[get1Dfrom3D(obase + 2*box - 1, 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &shftfp[get1Dfrom3D(lev-1-1, 0, 0, n-2, terms, terms)], phirs);
            
            for(int sc = 1; sc<=npm1;  sc++) {
                for(int term = 1; term<=p; term++) {
                    phi[get1Dfrom3D(base + box - 1, sc-1, term-1, 2*t + log2np - 3, p-1, terms)] = phils[get1Dfrom2D(sc-1, term-1, p-1, terms)] + phirs[get1Dfrom2D(sc-1, term-1, p-1, terms)];
                }
            }
        }
    }

    //At higher levels, communication is required
    for(int lev = log2np - 1; lev>=2; lev--) {
        int obase = base;
        base = base-1;
        int intrvl = p/pow(2, lev);
        int intvlo = intrvl/2;

        if(myid%intrvl == 0) {
            MPI_Recv(phr, pnpm1, MPI_DOUBLE_COMPLEX, nextd[get1Dfrom2D(lev+1-1, 0, p, 3)], lev, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            mpp(&phi[get1Dfrom3D(obase + 1 - 1 , 0, 0, 2*t + log2np - 3, p-1, terms)], terms, npm1, &shftfp[get1Dfrom3D(lev-1-1, 0, 0, n-2, terms, terms)], phils);
            mpp(phr, terms, npm1, &shftfp[get1Dfrom3D(lev-1-1, 0, 0, n-2, terms, terms)], phirs);

            for(int sc=1; sc<=npm1; sc++) {
                for(int term=1; term<=terms; terms++) {
                   phi[get1Dfrom3D(base + 1 - 1, sc-1, term-1, 2*t + log2np - 3, p-1, terms)] = phils[get1Dfrom2D(sc-1, term-1, p-1, terms)] + phirs[get1Dfrom2D(sc-1, term-1, p-1, terms)]; 
                }
            }
        }

        else if(myid % intrvl == intvlo) {
            MPI_Send(phi[get1Dfrom3D(obase + 1 - 1 , 0, 0, 2*t + log2np - 3, p-1, terms)], pnpm1, MPI_DOUBLE_COMPLEX, prevd[get1Dfrom2D(lev+1-1, 0, p, 3)], lev, MPI_COMM_WORLD);
        }
    }

    int tag = log2np;

    /* Pack:  phi, qr -> packns, packps */
    int pk = 0;

    for(int lev = lup+1; lev<=n; lev++) {
        int nl = pow(2, lev)/p;
        base = nl + log2np - 3;

        for(int sc=1; sc<=npm1; sc++) {
            for(int term=1; term<=terms; term++) {
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

    for(int sc=1; sc<=npm1; sc++) {
        for(int j=1; j<=b; j++) {
            pk++;
            packns[pk-1] = qr[get1Dfrom3D(t-1, sc+1-1, j-1, t, p, b)];
            packps[pk-1] = qr[get1Dfrom3D(0, sc+1-1, j-1, t, p, b)];
        }
    }

    /* Now pk .eq. (n - lup)*npm1*2*p + npm1*nieach */
    int lpkp = pk;

    if(t==1) {
        for(int sc = p/2; sc<=npm1; sc++) {
            int sce = sc - p/2 + 1;
            packnn[sce-1] = dotp(b, &qr[get1Dfrom3D(0, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 0, p/2, 3*b)]);
        }
    }
    else {
        for(int sc=p/2; sc<=npm1; sc++) {
            int sce = sc - p/2 + 1;
            pk++;
            packns[pk-1] = dotp(b, &qr[get1Dfrom3D(t-1-1, sc+1-1, 0, t, p, b)], &cotsh[get1Dfrom2D(sce-1, 0, p/2, 3*b)]);
        }
    }

    int lpkn = pk;

    /* Data transfer */
    tag++;
    MPI_Sendrecv(packps, lpkp, MPI_DOUBLE_COMPLEX, previd, tag, packnr, lpkp, MPI_DOUBLE_COMPLEX, nextid, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    tag++;
    MPI_Sendrecv(packns, lpkn, MPI_DOUBLE_COMPLEX, nextid, tag, packpr, lpkn, MPI_DOUBLE_COMPLEX, previd, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if(t==1) {
        tag++;
        MPI_Sendrecv(packnn, p/2, MPI_DOUBLE_COMPLEX, nextd[get1Dfrom2D(log2np, 1, p, 3)], tag, qcpp, p/2, MPI_DOUBLE_COMPLEX, prevd[get1Dfrom2D(log2np, 1, p, 3)], tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    /* Unpack:  packnr -> phin, qrn
          packpr -> phip, qrp, qcpp */

    int pk = 0;
    int lr = 0;

    for(int lev = lup+1; lev<=n; lev++) {
        lr++;
        nl = pow(2, lev)/p;
        base = nl + log2np - 3;

        for(int sc=1; sc<=npm1; sc++) {
            for(int term=1; term<=p; term++) {
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

    for(int sc=1; sc<=npm1; sc++) {
        for(int j=1; j<=b; j++) {
            pk++;
            qrn[get1Dfrom2D(sc-1, j-1, p, b)] = packnr[pk-1];
            qrp[get1Dfrom2D(sc-1, j-1, p, b)] = packpr[pk-1];
        }
    }

    /* Now pk .eq. (n - lup)*npm1*2*p + nieach*npm1 */
    if(t>=1) {
        for(int sc=1; sc<=p/2; sc++) {
            pk++;
            qcpp[sc-1] = packpr[pk-1];
        }
    }


}

void mft(int lq, double complex* qr, int dir, int p, int myid, int terms, int b, double* w, double* v, int sz1, int szwk) {
    MPI_Barrier(MPI_COMM_WORLD);
    int s = lq/b;
    int t = s/p;
    int n = my_log2(s);
    int log2np = my_log2(p);

    int ia[57];
    int ind = 0;

    for(int j=0; j<=7; j++) {
        ia[j] = ind;
        ind += terms*terms*(n-2);
    }
    ia[8] = ind;
    ind += p*p;
    ia[9] = ind;
    ind += terms*b;
    ia[10] = ind;
    ind += terms*(p - 1)*b;
    ia[11] = ind;
    ind += (terms*p)/2;
    ia[12] = ind;
    ind += (3*b*p)/2;
    ia[13] = ind;
    ind += 3*b*(b-1)*(p-1);
    ia[14] = ind;
    ind += (3*b*p)/2;
    ia[15] = ind;
    ind += terms;
    ia[16] = ind;
    ind += 2*(p-1)*terms;

    ind = 0;
    for(int j=17; j<=18; j++) {
        ia[j] = ind;
        ind+= 2*terms*(p-1)*(2*t + log2np - 3);
    }
    for(int j=21; j<=22; j++) {
        ia[j] = ind;
        ind += 2*p;
    }
    for(int j=23; j<=25; j++) {
        ia[j] = ind;
        ind += 2*terms*(p-1);
    }
    for(int j=26; j<=27; j++) {
        ia[j] = ind;
        ind += 2*(p-1)*(b + n*2*terms);
    }
    for(int j=28; j<=29; j++) {
        ia[j] = ind;
        ind += 2*((p-1)*(b + n*2*terms) + p/2);
    }
    ia[30] = ind;
    ind += 2*b*p;
    ia[31] = ind;
    ind += 2*b*(p-1);
    for(int j=32; j<=33; j++) {
        ia[j] = ind;
        ind+= 2 * (p/2);
    }
    for(int j=34; j<=43; j++) {
        ia[j] = ind;
        ind += 2*terms*(p-1);
    }
    for(int j=44; j<=45; j++) {
        ia[j] = ind;
        ind += 2 * terms * (p-1) * 2 * n;
    }
    ia[46] = ind;
    ind += 2*lq;
    ia[47] = ind;
    ind += 2*(p-1);
    for(int j=48; j<=49; j++) {
        ia[j] = ind;
        ind += 2 * (p/2);
    }
    ia[50] = ind;
    ind = ind + b;
    for(int j=51; j<=53; j++) {
        ia[j] = ind;
        ind += terms;
    }
    for(int j=54; j<=55; j++) {
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
        &v[ia[49]], &v[ia[50]], &v[ia[51]], &v[ia[52]], &v[ia[53]], &v[ia[54]], &v[ia[55]]
    );
}