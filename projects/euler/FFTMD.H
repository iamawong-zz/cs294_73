#ifndef _FFTMD_H_
#define _FFTMD_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <list>
#include <tr1/memory>
#include "Box.H"
#include "MDArray.H"
#include "FFT1D.H"
using namespace std;
class FFTMD
{
public:
  FFTMD();
  FFTMD(std::tr1::shared_ptr<FFT1D> a_fft1dptr);
  void forwardCC(MDArray<complex<double> > & a_f) const;
  void inverseCC(MDArray<complex<double> > & a_fHat) const;
  const int& getN() const;
  const int& getM() const;
private:
  int m_N;
  int m_M;
  std::tr1::shared_ptr<FFT1D> m_fft1dPtr;
};
#endif
