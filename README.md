# fmm-fft

Parallel FMM Accellerated FFT on CPUs using OpenMPI.
Based on the works of [Edelman et. al.](http://www-math.mit.edu/%7Eedelman/homepage/papers/edelman97future.pdf).

### Authors
* [Amarnath Karthi](https://github.com/97amarnathk)
* [Chahak Mehta](https://github.com/chahak13)

## Introduction

Traditional distributed Fast Fourier Transforms (FFTs) require 3 _all-to-all_ communications which cause severe communication bottlenecks.
FMM-FFT uses the Fast Multipole Method to approximate the solutions of the FFT, thereby reducing the communication required.

## Dependencies

* FFTW3
* OpenMPI
* OpenMP

## Codes

* [src/fftw_mpi.c](https://github.com/97amarnathk/fmm-fft/blob/master/src/fftw_mpi.c)  : Distributed FFT
* [src/fmmfft.c](https://github.com/97amarnathk/fmm-fft/blob/master/src/fmmfft.c)    : Distributed FMM Accellerated FFT
* [src/fmmfft_error.c](https://github.com/97amarnathk/fmm-fft/blob/master/src/fmmfft_error.c)  : Comparison of FFT and FMM-FFT outputs

## Installation
1. Edit `makefile` and change the `-L` and `-I` in the `local` target.
2. `make local`
3. Run the out file.

## References
[1] [The Future Fast Fourier Transform?](http://www-math.mit.edu/%7Eedelman/homepage/papers/edelman97future.pdf), by Edelman et. al. , SIAM Journal on Scientific Computing. SIAM J Scientific Computing 20 1094-1114 (1998)

[2] [Re-Introduction of Communication-AvoidingFMM-Accelerated FFTs with GPU Acceleration](http://www.harperlangston.com/papers/hpec13.pdf), by Langston et. al. , IEEE Conference on High Performance Extreme Computing (HPEC '13)


### TODO:
* Add headers for FMM-FFT functions.