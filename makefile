local : 
	mpicc src/fmmfft.c -lfftw3 -lm -w -o fmmfft.out -O2
	mpicc src/fftw_serial.c -lfftw3 -lm -w -o fftw_serial.out 
	mpicc src/fftw_mpi.c -lm -lfftw3 -lfftw3_mpi -o fft_mpi.out -O2

cluster :
	mpicc src/fmmfft.c -lfftw3 -lm -w -L/home/201501005/fftw3l/lib/ -I/home/201501005/fftw3l/include/ -std=c99 -o fmmfft.out -O2
	mpicc src/fftw_serial.c -lfftw3 -lm -w -L/home/201501005/fftw3l/lib/ -I/home/201501005/fftw3l/include/ -std=c99 -o fftw_serial.out -O2
	mpicc fftw_mpi.c -lm -lfftw3 -L/home/201501005/fftw3l/lib/ -I/home/201501005/fftw3l/include/ -lfftw3_mpi -O2 -std=c99