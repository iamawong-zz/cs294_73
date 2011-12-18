#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

system := $(shell uname)
FFTWDIR = $(HOME)/fftw
FFTINCFLAGS = -I$(FFTWDIR)/include
VPATH = ../../resources/VisitWriter
CFLAGS = -g -Wall
#CFLAGS = -O3  
LIBS:=  -L$(FFTWDIR)/lib -lfftw3 
CXX = g++ 

CFLAGS+=$(FFTINCFLAGS) -I../../resources/VisitWriter

MOBJS= PowerItoI.o FFTMDSineTest.o FFTW1DSine.o FFTMDSine.o Box.o WriteMDArray.o VisitWriter.o ParticleSet.o

POISOBJS= PowerItoI.o FFTW1DSine.o FFTMDSine.o Box.o WriteMDArray.o VisitWriter.o PoissonInfinite.o PoissonTest.o

POBJS= ParticleTest.o ParticleSet.o Box.o WriteMDArray.o VisitWriter.o

VOBJS = PowerItoI.o FFTW1DSine.o FFTMDSine.o MLCVortexMethod.o ParticleSet.o ParticleVelocities.o  PoissonInfinite.o Box.o VortexTest.o WriteMDArray.o VisitWriter.o 

ParticleTest: $(POBJS)
	$(CXX) $(CFLAGS) -o ParticleTest $(POBJS)


testvortex: $(VOBJS) GNUmakefile
	$(CXX) $(CFLAGS) -o testvortex $(VOBJS)  $(LIBS)


%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

clean:
	rm *.o *.d

cleanvtk: 
	rm *.vtk

# dependency info for *existing* .o files
-include $(POBJS:.o=.d)
-include $(VOBJS:.o=.d)
