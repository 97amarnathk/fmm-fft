local : 
	mpicc src/converted.c -lfftw3 -lm -w -o fmmfft.out

cluster :
	mpicc src/converted.c -lfftw3 -lm -w -L/home/201501005/fftw3l/lib/ -I/home/201501005/fftw3l/include/ -std=c99 -o fmmfft.out