#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "MDArray.H"
#include "FFT1D.H"
#include "FFTMD.H"
using namespace std;
FFTMD::FFTMD(FFT1D* a_fft1dPtr)
{
  m_fft1dPtr = a_fft1dPtr;
  m_M = m_fft1dPtr->getM();
  m_N = m_fft1dPtr->getN();
} 
void FFTMD::forwardCC(MDArray<complex<double> > & a_f) const
{
  int low[DIM],high[DIM],tuple[DIM];
  vector<complex<double> > f1d(m_N);
  vector<complex<double> > fHat1d(m_N);
  
  for (int dir = 0;dir < DIM ; dir++)
    {
      for (int dir2 = 0;dir2 < DIM;dir2++)
        {
          low[dir2]= 0;
          high[dir2] = m_N-1;
        }
      high[dir]=0;
      Box base(low,high);
      for (int k = 0;k < base.sizeOf();k++)
        {
          base.tupleIndex(k,tuple);
          for (int l = 0 ; l < m_N;l++)
            {
              tuple[dir]= l;
              f1d[l] = a_f[tuple];
            }
          m_fft1dPtr->forwardFFTCC(fHat1d,f1d);
          tuple[dir] = 0;
          for (int l = 0 ; l < m_N;l++)
            {
              tuple[dir] = l;
              a_f[tuple] = fHat1d[l];
            }
        }
    }     
};
void FFTMD::inverseCC(MDArray<complex<double> > & a_fHat) const
{int low[DIM],high[DIM],tuple[DIM];
  vector<complex<double> > f1d(m_N);
  vector<complex<double> > fHat1d(m_N);
  
  for (int dir = 0;dir < DIM ; dir++)
    {
      for (int dir2 = 0;dir2 < DIM;dir2++)
        {
          low[dir2]= 0;
          high[dir2] = m_N-1;
        }
      high[dir]=0;
      Box base(low,high);
      for (int k = 0;k < base.sizeOf();k++)
        {
          base.tupleIndex(k,tuple);
          for (int l = 0 ; l < m_N;l++)
            {
              tuple[dir] = l;
              fHat1d[l] = a_fHat[tuple];
            }
          m_fft1dPtr->inverseFFTCC(f1d,fHat1d);
          tuple[dir] = 0;
          for (int l = 0 ; l < m_N;l++)
            {
              tuple[dir] = l;
              a_fHat[tuple] = f1d[l];
            }
        }
    }     
};
const int& FFTMD::getN() const
{
  return m_N;
  
};
const int& FFTMD::getM() const
{
  return m_M;
};
