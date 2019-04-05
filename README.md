# fmm-fft
Parallel FMM Accellerated FFT on CPUs


## Functions

### Verified

### Written

* mfti
* mftii
* initg
* flip
* shftf
* shftl
* vx
* wx
* evlmf
* mft

### TODO

* mftint


## Verify the following

* All equations. For example `a / b / c` is correctly written.

### Arrays
* All computations are done in 1 base indexing. All arrays are stored in 0 based indexing. For example:

```C
// All arrays should be used like this
for(int j=1; j<=terms; j++) {
  double th = f(j);
  arr[j-1] = temp;
}
```

* Fortran is column major. ie `A(2, 1, 10)` is placed next to `A(1, 1, 10)`. In C, `A(1, 2, 10)` is the neighbor of `A(1, 1, 10)`. Therefore if a fortran array ARR has dimensions I x J x K ,  it will be stored in C as `A[K][J][I]` array. 

As an example consider the following fortran to C conversion. Fortran code is:

```fortran
double precision ARR(I, J, K)

do i=1, I
  do j=1, J
    do k=1, K
      ARR(i, j, k) = f(i, j, k)
     end do
  end do
end do
```

It should be converted to the following C code:

```C
// stored as transpose, to emulate column major indexing
double ARR[K][J][I];

// loops in 1 base indexing
for(int i=1; i<=I; i++) {
  for(int j=1; j<=J; j++) {
    for(int k=1; k<=K; k++) {
      // stored in 0 based indexing (subtrct 1), store as transpose.
      ARR[k-1][j-1][i-1] = f(i, j, k); 
    }
  }
}
```
