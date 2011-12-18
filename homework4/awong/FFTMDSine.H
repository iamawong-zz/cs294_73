#ifndef _FFTMDSINE_H_
#define _FFTMDSINE_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "PowerItoI.H"
#include "MDArray.H"
#include "FFT1D.H"
//#include "FFTW1D.H"
#include "FFTW1DSine.H"
using namespace std;
class FFTMDSine
{
public:
  FFTMDSine(unsigned int a_M);
  ~FFTMDSine();
  void transform(MDArray<double> & a_f) const;
  const int& getN() const;
private:
  int m_N;
  FFTW1DSine* m_fft1d;
  Box m_box;
};
#endif
