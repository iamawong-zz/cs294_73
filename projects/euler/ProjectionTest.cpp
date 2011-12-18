#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "MDArray.H"
#include "FFT1D.H"
#include "FFTW1D.H"
#include "FFTMD.H"
#include "PoissonSolver.H"
#include "Projection.H"
#include "FieldData.H"
#include "DeltaVelocity.H"
#include "RK4.H"
#include "WriteMDArray.H"
#include "PowerItoI.H"

int main(int argc, char* argv[])
{
  // Make N dependent on M (for safety)
  unsigned int M = 7;
  unsigned int N = Power(2, M);

  std::tr1::shared_ptr<FFT1D> fft1dptr(new FFTW1D(M));
  FieldData velocities(fft1dptr, 1);
  Projection testing(fft1dptr);

  int low[DIM] = {0,0};
  int high[DIM] = {N-1,N-1};
  Box bx(low,high);
  double rho = 1.0/30.0; // careful when doing floating point arithmetic. don't want truncated int values
  double delta = 0.05;
  double dx = 1.0/((double)N); // careful when doing floating point arithmetic. don't want truncated int values

  for(int i=0; i<bx.sizeOf(); i++) {
    int tuple[DIM];
    bx.tupleIndex(i, tuple);

    if(tuple[1]*dx <= 0.5) {
      velocities.m_data[0][tuple] = tanh((tuple[1]*dx-0.25)/rho);
    } else if(tuple[1]*dx > 0.5) {
      velocities.m_data[0][tuple] = tanh((0.75-tuple[1]*dx)/rho);
    }

    velocities.m_data[1][tuple] = delta*sin(2.0*M_PI*tuple[0]*dx);
  }

  testing.applyProjection(velocities);
  MDArray<Real> divergenceSc(bx); // need to set the box for the MDArray
  testing.divergence(divergenceSc, velocities);

  MDWrite("ProjectionTest",divergenceSc);

  ///somehow print out velocities
}

// Stuff from Projection:
//  Projection(std::tr1::shared_ptr<FFT1D> a_fft1dPtr);
//  void applyProjection(FieldData& a_velocity);

