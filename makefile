local : 
	mpicc src/fmmfft.c -lfftw3 -lm -w -o fmmfft.out -O2

cluster :
	mpicc src/fmmfft.c -lfftw3 -lm -w -L/home/201501005/fftw3l/lib/ -I/home/201501005/fftw3l/include/ -std=c99 -o fmmfft.out -O2