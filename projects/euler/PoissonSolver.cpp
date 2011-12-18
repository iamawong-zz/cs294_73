#include <iostream>
#include "PoissonSolver.H"
using namespace std;

PoissonSolver::PoissonSolver()
{
}
PoissonSolver::PoissonSolver(std::tr1::shared_ptr<FFT1D> a_fft1dPtr)
{
  m_fft1dptr = a_fft1dPtr;
  m_N = a_fft1dPtr->getN();
  m_M = a_fft1dPtr->getM();
}

void PoissonSolver::solve(MDArray<Real>& a_Rhs) const {
  //create a contructor
  Box bx = a_Rhs.getBox(); 
  MDArray<complex<double> > fftwForward(bx);
  MDArray<complex<double> > fourierCoef(bx);

  double h = 1.0/((double)m_N);
   h = .0001;
  //double h = 1/m_N; // integer division, becomes 0.0 in the end!

  for (int k = 0; k < bx.sizeOf();k++) {
    fftwForward[k] = complex<double>(a_Rhs[k], 0.0);
  }
  FFTMD fftmd = FFTMD(m_fft1dptr);
  fftmd.forwardCC(fftwForward);

  for (int k=0; k<bx.sizeOf(); k++) {
    int index[2];
    bx.tupleIndex(k,index);
    int i = index[0];
    int j = index[1];
    if (k == 0) {
      fourierCoef[index] = complex<double>(0.0,0.0);
    } else {
      complex<double> div((2.0*cos(2.0*M_PI*((double)i)*h) - 2.0*cos(2.0*M_PI*((double)j)*h) - 4.0)/h/h,0.0);
      fourierCoef[index] =fftwForward[index]/div;
    }
  }

  fftmd.inverseCC(fourierCoef);

  for (int k=0; k<bx.sizeOf(); k++) {
    a_Rhs[k] = real(fourierCoef[k]) * h * h;
  }

}
