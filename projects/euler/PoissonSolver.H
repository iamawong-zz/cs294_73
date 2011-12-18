#ifndef _POISSONSOLVER_H_
#define _POISSONSOLVER_H_

#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <complex>
#include "Box.H"
#include "MDArray.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include "Real.H"
using namespace std;

class PoissonSolver
{
public:
  PoissonSolver();
  PoissonSolver(std::tr1::shared_ptr<FFT1D> a_fft1dPtr);
  void solve(MDArray<Real>& a_Rhs) const;
  int m_M,m_N;
  std::tr1::shared_ptr<FFT1D> m_fft1dptr;
};

#endif

