#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
using namespace std;
#include "PowerItoI.H"
#include "FFT1DBRI.H"
#include "FFT1DRecursive.H"
#include "FFTW1D.H"
#include "MDArray.H"
#include "Box.H"
#include "FFTMD.H"
double test(const FFTMD& a_fftmd)
{
  int N = a_fftmd.getN();
  int low[DIM],high[DIM],tuple[DIM];
  
  for (int dir = 0;dir <DIM;dir++)
    {
      low[dir] = 0;
      high[dir] = N-1;
    }
  Box b(low,high);
     
  MDArray<complex<double> > f(b);
  MDArray<complex<double> > fSave(b);
  MDArray<complex<double> > fHat(b);
  double h = 1./N;
  
  for (int k=0; k < b.sizeOf();k++)
    {
      b.tupleIndex(k,tuple);
      double x = 0.;
      
      for (int dir = 0;dir < DIM;dir++)
        {
          x += pow(tuple[dir]*h,2);
        }
      x *=32*32;
      f[k] = complex<double >(exp(-x),0);
      //f[k] = complex<double>(1.,0.);
      fSave[k] = f[k];
    }
  a_fftmd.forwardCC(f);
  a_fftmd.inverseCC(f);
  double maxerr = 0.;
    int normalize = Power(N,DIM);
  
  for (int j =0; j < b.sizeOf() ; j++)
      {
        if (fabs(real(f[j])/normalize - real(fSave[j])) > maxerr )
          {
            maxerr = fabs(real(f[j])/normalize - real(fSave[j]));
          }
      }
  return maxerr;
}
int main(int argc, char* argv[])
{
  int M;
  // int inputMode;
  
  sscanf(argv[1],"%d",&M);
  // sscanf(argv[2],"%d",&inputMode);
  //FFT1DRecursive fft1d(M);	
  //FFT1DBRI fft1d(M);
  FFTW1D fft1d(M);
  
  FFT1D* p_fft = dynamic_cast<FFT1D*>(&fft1d);
  FFTMD fftmd(p_fft);
  double error = test(fftmd);
  cout << "test 1: error in Gaussian  = " << error << endl;
};
