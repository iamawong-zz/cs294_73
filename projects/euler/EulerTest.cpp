#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "MDArray.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include "FFTW1D.H"
#include "PoissonSolver.H"
#include "Projection.H"
#include "FieldData.H"
#include "DeltaVelocity.H"
#include "RK4.H"
#include "WriteMDArray.H" 
#include "ComputeEulerRHS.H"

double returnMax(MDArray<double> velocities[DIM]) {
  double umax = 0.;
  double vmax = 0.;
  Box grids = velocities[0].m_box;
  
  for(int i=0; i<grids.sizeOf(); i++) {
    int tuple[DIM];
    double currVelocities[DIM];
    grids.tupleIndex(i, tuple);
    currVelocities[0] = velocities[0][tuple];
    currVelocities[1] = velocities[1][tuple];
    umax = max(currVelocities[0], umax);
    vmax = max(currVelocities[1], vmax);
  }

  return max(umax, vmax);
} 

void calcVorticity(MDArray<double> &vorticity, const Box bx, FieldData &velocities, double dx) {
    velocities.fillGhosts();
  for(int i=0; i<bx.sizeOf(); i++) {
    int tuple[DIM];
    int tupleShiftLow[DIM];
    int tupleShiftHigh[DIM];
    double dudy;
    double dvdx;
    bx.tupleIndex(i,tuple);

    tupleShiftLow[0] = tuple[0];
    tupleShiftLow[1] = tuple[1]-1;
    tupleShiftHigh[0] = tuple[0];
    tupleShiftHigh[1] = tuple[1]+1;
    dudy = (velocities[0][tupleShiftHigh] - velocities[0][tupleShiftLow])/(2*dx);
    tupleShiftLow[0] = tuple[0]-1;
    tupleShiftLow[1] = tuple[1];
    tupleShiftHigh[0] = tuple[0]+1;
    tupleShiftHigh[1] = tuple[1];	
    dvdx = (velocities[1][tupleShiftHigh] - velocities[1][tupleShiftLow])/(2*dx);
    vorticity[tuple] = dudy-dvdx;
  }
}

int main(int argc, char* argv[])
{
  // Make N dependent on M (for safety)
  unsigned int M = 7;
  unsigned int N = Power(2, M);

  std::tr1::shared_ptr<FFT1D> fft1dptr(new FFTW1D(M));
  FieldData velocities(fft1dptr, 1);
  int low[DIM] = {0,0};
  int high[DIM] = {N-1,N-1};
  Box bx(low,high);
  double rho = 1.0/30.0; // careful when doing floating point arithmetic. don't want truncated int values
  double delta = 0.05;
  double dx = 1.0/((double)N); // careful when doing floating point arithmetic. don't want truncated int values
  double timeEnd = 0.8;
  double C = 0.5;
  char name[10];
  int nstop = 15;

  /// Initial Conditions
  for(int i=0; i<bx.sizeOf(); i++) {
    int tuple[DIM];
    bx.tupleIndex(i, tuple);
    velocities.m_data[1][tuple] = delta*sin(2.0*M_PI*tuple[0]*dx);
    if(tuple[1]*dx <= 0.5) {
      velocities.m_data[0][tuple] = tanh((tuple[1]*dx-0.25)/rho);
    } else if(tuple[1]*dx > 0.5) {
      velocities.m_data[0][tuple] = tanh((0.75-tuple[1]*dx)/rho);
    }
  }

  RK4<FieldData, ComputeEulerRHS, DeltaVelocity> integrator;

  double dt = C*dx/returnMax(velocities.m_data);
  double time=0.0;

      MDArray<double> vorticity(bx);
      calcVorticity(vorticity, bx, velocities, dx);
      sprintf(name, "Vorticity.%.4f.%.4f", timeEnd, time);
      MDWrite(name, vorticity);
      int nstep = 0;
      while(nstep < nstop)
    {
      integrator.advance(time, dt, velocities);
      time += dt;
      nstep++;
      calcVorticity(vorticity, bx, velocities, dx);
      sprintf(name, "Vorticity.%d",nstep);
      MDWrite(name, vorticity);
      cout << "dt = " << dt << endl;
      dt = C*dx/returnMax(velocities.m_data);
    }
}
