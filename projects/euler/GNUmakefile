#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

system := $(shell uname)
FFTWDIR = $(HOME)/fftw
FFTINCFLAGS = -I$(FFTWDIR)/include
#for Johnny and possibly Ivan (Ubuntu machines)
ifeq ((system),Linux)
  #LAPACKFLAGS = -lcblas -lclapack -I/usr/include/atlas -L/usr/local/lib -L/usr/lib
  LAPACKFLAGS = -lcblas -llapack -I/usr/include/atlas
endif
# for OSX
ifeq ($(system),Darwin)
  LAPACKFLAGS = -framework Accelerate
endif
# don't actually need LAPACK
LAPACKFLAGS = 

VPATH = ../../resources/VisitWriter
CFLAGS = -g -Wall $(LAPACKFLAGS)
#CFLAGS = -O3
LIBS:=  -L$(FFTWDIR)/lib -lfftw3 
CXX = g++ 

CFLAGS+=$(FFTINCFLAGS) -I../../resources/VisitWriter

REQOBJS= FFTMD.o FFTW1D.o PowerItoI.o Box.o WriteMDArray.o PoissonSolver.o VisitWriter.o

EULEROBJS= Projection.o AdvectionOperator.o ComputeEulerRHS.o

ProjectionTestObjs = $(REQOBJS) FieldData.o DeltaVelocity.o Projection.o

AdvectionTestObjs = $(REQOBJS) FieldData.o DeltaVelocity.o Projection.o AdvectionOperator.o

ComputeEulerRHSTestObjs = $(REQOBJS) FieldData.o DeltaVelocity.o Projection.o AdvectionOperator.o ComputeEulerRHS.o

FFTMD.o: FFTMD.cpp FFTMD.H Box.H MDArray.H MDArrayImplem.H GNUmakefile
	$(CXX) $(CFLAGS) -c FFTMD.cpp

FFTW1D.o: FFTW1D.H FFTW1D.cpp GNUmakefile
	$(CXX) $(CFLAGS) $(FFTINCFLAGS) -c -o FFTW1D.o FFTW1D.cpp

PowerItoI.o: PowerItoI.cpp PowerItoI.H  GNUmakefile
	$(CXX) $(CFLAGS) -c PowerItoI.cpp

Box.o: Box.cpp Box.H  GNUmakefile
	$(CXX) $(CFLAGS) -c Box.cpp

WriteMDArray.o: WriteMDArray.cpp WriteMDArray.H MDArray.H MDArrayImplem.H VisitWriter.H GNUmakefile
	$(CXX) $(CFLAGS) -c WriteMDArray.cpp

DeltaVelocity.o: DeltaVelocity.H DeltaVelocity.cpp FieldData.H Box.H MDArray.H MDArrayImplem.H GNUmakefile
	$(CXX) $(CFLAGS) -c DeltaVelocity.cpp

FieldData.o: FieldData.H FieldData.cpp DeltaVelocity.H Box.H MDArray.H MDArrayImplem.H FFT1D.H GNUmakefile
	$(CXX) $(CFLAGS) -c FieldData.cpp

PoissonSolver.o: PoissonSolver.H PoissonSolver.cpp Box.H FFT1D.H FFTMD.H GNUmakefile
	$(CXX) $(CFLAGS) -c PoissonSolver.cpp

Projection.o: Projection.H Projection.cpp PoissonSolver.H Box.H MDArray.H MDArrayImplem.H FieldData.H DeltaVelocity.H GNUmakefile
	$(CXX) $(CFLAGS) -c Projection.cpp

AdvectionOperator.o: AdvectionOperator.H AdvectionOperator.cpp Box.H MDArray.H FieldData.H DeltaVelocity.H GNUmakefile
	$(CXX) $(CFLAGS) -c AdvectionOperator.cpp

MyMachineName:
	@echo $(system) $(HOME)

AdvectionTest: $(AdvectionTestObjs) AdvectionTest.cpp GNUmakefile
	$(CXX) $(CFLAGS) AdvectionTest.cpp $(AdvectionTestObjs) -o AdvectionTest $(LIBS)

ProjectionTest: $(ProjectionTestObjs) ProjectionTest.cpp GNUmakefile
	$(CXX) $(CFLAGS) ProjectionTest.cpp $(ProjectionTestObjs) -o ProjectionTest $(LIBS)

EulerTest: $(ComputeEulerRHSTestObjs) EulerTest.cpp RK4.H GNUmakefile
	$(CXX) $(CFLAGS) EulerTest.cpp $(ComputeEulerRHSTestObjs) -o EulerTest $(LIBS)

%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

clean:
	rm *.o *.d
cleanExecs:
	rm AdvectionTest ProjectionTest EulerTest
