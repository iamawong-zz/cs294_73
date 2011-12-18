#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "MDArray.H"
#include "FFTW1DSine.H"
#include "FFTMDSine.H"
using namespace std;
FFTMDSine::FFTMDSine(unsigned int a_N)
{
  m_fft1d = new FFTW1DSine(a_N);
  m_N = a_N;
  int lo[DIM] ={0,0};
  int hi[DIM] ={m_N+1,m_N+1};
  m_box=Box(lo,hi);
}
void FFTMDSine::transform(MDArray<double> & a_f) const
{
  int low[DIM],high[DIM],tuple[DIM];
  vector<double>   f1d(m_N);
  vector<double>  fHat1d(m_N);
  int lowCorner[DIM];
  a_f.getBox().getLowCorner(lowCorner);
  for (int dir = 0;dir < DIM; dir++)
    {
      a_f.shift(dir,-lowCorner[dir]);
    }
  
  
  for (int dir = 0;dir < DIM ; dir++)
    {
      for (int dir2 = 0;dir2 < DIM;dir2++)
        {
          low[dir2]= 1;
          high[dir2] = m_N+1;
        }
      high[dir]=1;
      Box base(low,high);
      for (int k = 0;k < base.sizeOf();k++)
        {
          base.tupleIndex(k,tuple);
          for (int l = 0 ; l < m_N;l++)
            {
              tuple[dir]= l+1;
              f1d[l] = a_f[tuple];
            }
          m_fft1d -> transform(fHat1d,f1d);
          for (int l = 0 ; l < m_N;l++)
            {
              tuple[dir] = l+1;
              a_f[tuple] = fHat1d[l]/sqrt(2*(m_N+1));
            }
        }
    }
  for (int dir = 0; dir < DIM;dir++)
    {
      a_f.shift(dir,lowCorner[dir]);
    }
};
const int& FFTMDSine::getN() const
{
  return m_N;
};
FFTMDSine::~FFTMDSine(){delete m_fft1d;};
