CC  = gcc -Wall -m64 -march=k8 -msse3 -mfpmath=sse -O3
#F90 = ifort -warn -nofor_main -ipo -O3 -no-prec-div -static -openmp -assume cc_omp -complex_limited_range \
            -debug extended -traceback
F90 = gfortran -m64 -mavx -O3 -march=native -mtune=native -fopenmp -floop-block \
      -ffast-math -fcx-fortran-rules -mrecip -fblas-matmul-limit=50 \
      -fno-realloc-lhs -fbacktrace -g -ffree-line-length-none
F90L = $(F90)
#LAPACK = -L /usr/local/intel/icte/3.1.1/mkl/10.0.3.020/lib/em64t -lmkl_lapack -lmkl_em64t -lguide -lpthread
#LAPACK = -L /usr/local/intel/ics/2013.0.028/mkl/lib/intel64 -lmkl_lapack -lmkl_em64t -lguide -lpthread 
LAPACK =  -L/usr/lib64 -lblas -llapack 
FFTW3 =  -L/usr/lib64 -lfftw3 -lfftw3f
#DXINC = -DUSE_DX -I/home/software/OpenDX-4.4.4/dx/include
#DXLIB = -L/home/software/OpenDX-4.4.4/dx/lib_linux/ -lDXlite
DXINC =
DXLIB =

fftw.o:	fftw_fftw3.f90 accuracy.o
	$(F90) -c fftw_fftw3.f90 -o fftw.o
