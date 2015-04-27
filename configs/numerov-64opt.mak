CC  = gcc -Wall -m64 -march=k8 -msse3 -mfpmath=sse -O3
F90 = ifort -warn -nofor_main -ipo -O3 -no-prec-div -static -xP -openmp -assume cc_omp -complex_limited_range \
            -debug extended -traceback 
F90L = $(F90)
LAPACK = -L /usr/local/intel/icte/3.1.1/mkl/10.0.3.020/lib/em64t -lmkl_lapack -lmkl_em64t -lguide -lpthread
FFTW3 = -L /usr/local/fftw/3.2/intel/lib/ -lfftw3 -lfftw3f
DXINC = -DUSE_DX -I/home/software/OpenDX-4.4.4/dx/include
DXLIB = -L/home/software/OpenDX-4.4.4/dx/lib_linux/ -lDXlite

fftw.o:	fftw_fftw3.f90 accuracy.o
	$(F90) -c fftw_fftw3.f90 -o fftw.o
