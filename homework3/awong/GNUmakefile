#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

system := $(shell uname)
FFTWDIR = $(HOME)/fftw
FFTINCFLAGS = -I$(FFTWDIR)/include
CFLAGS = -g -Wall
#CFLAGS = -O3
BLASDIR = $(HOME)/OpenBLAS
BLASINCFLAGS = -I$(BLASDIR)/include
BLASLIBFLAGS = -L$(BLASDIR)/lib -lopenblas -lgfortran
ifeq ($(system),Darwin)
  BLASFLAGS = -framework Accelerate
endif
#CFLAGS = -O3  
LIBS:=  -L$(FFTWDIR)/lib -lfftw3 
CXX = g++ 


test1d: FFT1DBRI.cpp FFT1DRecursive.cpp FFT1DTest.cpp FFTCTBRI.cpp PowerItoI.cpp FFT1DBRI.o FFT1DTest.o PowerItoI.o FFTCTBRI.o FFT1DRecursive.o FFTW1D.o GNUmakefile
	  $(CXX) $(CFLAGS) $(FFTINCFLAGS) FFT1DBRI.o FFT1DTest.o PowerItoI.o FFTCTBRI.o FFT1DRecursive.o FFTW1D.o $(LIBS) -o test1d.ex

FFT1DTest.o: FFT1DTest.cpp FFT1D.H FFT1DBRI.cpp FFTCTBRI.H PowerItoI.H  GNUmakefile
	$(CXX) $(CFLAGS) $(FFTINCFLAGS) -c -o FFT1DTest.o FFT1DTest.cpp 

FFTW1D.o: FFTW1D.H FFTW1D.cpp GNUmakefile
	$(CXX) $(CFLAGS) $(FFTINCFLAGS) -c -o FFTW1D.o FFTW1D.cpp

FFTMDTest.o: FFTMDTest.cpp FFT1D.H FFTMD.H FFTMD.cpp FFT1DBRI.cpp FFTCTBRI.H PowerItoI.H MDArray.H MDArrayImplem.H  GNUmakefile
	$(CXX) $(CFLAGS) -c -o FFTMDTest.o FFTMDTest.cpp

FFTMD.o: FFTMDTest.cpp MDArray.H FFT1D.H FFTMD.H FFTMD.cpp FFT1DBRI.cpp FFTCTBRI.H PowerItoI.H MDArray.H  GNUmakefile
	$(CXX) $(CFLAGS) -c FFTMD.cpp

FFT1DBRI.o: FFT1D.H FFT1DBRI.cpp FFTCTBRI.H PowerItoI.H FFT1DBRI.H  GNUmakefile
	$(CXX) $(CFLAGS) -c FFT1DBRI.cpp 

FFT1DRecursive.o: FFT1D.H FFT1DRecursive.cpp FFT1DRecursive.H PowerItoI.H  GNUmakefile
	$(CXX) $(CFLAGS) -c -o FFT1DRecursive.o FFT1DRecursive.cpp 

FFTCTBRI.o: FFTCTBRI.cpp FFTCTBRI.H PowerItoI.H  GNUmakefile
	$(CXX) $(CFLAGS) -c FFTCTBRI.cpp

PowerItoI.o: PowerItoI.cpp PowerItoI.H  GNUmakefile
	$(CXX) $(CFLAGS) -c -o PowerItoI.o PowerItoI.cpp

Box.o: Box.cpp Box.H  GNUmakefile
	$(CXX) $(CFLAGS) -c Box.cpp

testMD: FFTMDTest.o FFT1DBRI.o FFTCTBRI.o PowerItoI.o Box.o FFTMD.o MDArrayImplem.H MDArray.H FFT1D.H GNUmakefile
	$(CXX) $(FFTINCFLAGS) FFTMDTest.o FFTW1D.o FFT1DRecursive.o FFT1DBRI.o FFTCTBRI.o PowerItoI.o Box.o FFTMD.o $(LIBS) -o testMD.ex

clean:
	rm *.o *.ex
naive: benchmark.cpp dgemm-naive.cpp GNUmakefile
	$(CXX) $(CFLAGS) $(BLASINCFLAGS) $(BLASFLAGS) -o naive benchmark.cpp dgemm-naive.cpp $(BLASLIBFLAGS)
blas: benchmark.cpp dgemm-blas.cpp GNUmakefile
	$(CXX) $(CFLAGS) $(BLASINCFLAGS) $(BLASFLAGS) -o blas benchmark.cpp dgemm-blas.cpp $(BLASLIBFLAGS)
