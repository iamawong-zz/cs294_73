#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
using namespace std;
#include "PowerItoI.H"
#include "FFT1D.H"
#include "FFTW1D.H"

#include "FFT1DRecursive.H"
#include "FFT1DBRI.H"
double test1(FFT1D* a_fftPtr)
{
  int N = a_fftPtr->getN();
  vector<complex<double > > f(N);
  vector<double > fSave(N);
  vector<complex<double> > fHat(N);
  double h = 1./N;
  
  for (int j=0; j<N;j++)
    {
      double x = (j*h - .5)*32;
      f[j] = complex<double>(exp(-pow(x,2)),0.);
      fSave[j] = real(f[j]);
    }
  a_fftPtr->forwardFFTCC(fHat,f);
  a_fftPtr->inverseFFTCC(f,fHat); 
  
  /*  for (int j = 0 ; j< N;j++)
    {
      int jout = j;
      if (jout >N/2) 
        {
          jout = j - N;
        }
       cout << j << " " << real(f[j])/N << " " << fSave[j] << endl;
       } */

  double maxerr = 0.;
  for (int j =0; j < N ; j++)
      {
        if (fabs(real(f[j])/N - fSave[j]) > maxerr )
          {
            maxerr = fabs(real(f[j])/N - fSave[j]);
          }
      }
  return maxerr;
}
int test2(FFT1D* a_fftPtr, int a_mode)
{
  int N = a_fftPtr->getN();
  vector<complex<double> > f(N);
  vector<complex<double> > fHat(N);
  double h = 1./N;
  complex<double> z(cos(2*M_PI*h*a_mode),sin(2*M_PI*h*a_mode));
  complex<double> zToTheJ(1.,0.);

  for (int j = 0;j < N;j++)
    {
      f[j] = zToTheJ;
      zToTheJ *= z;
    }
  
  a_fftPtr->forwardFFTCC(fHat,f);
  int numModes = 0;
  int maxMode = -1;
  
  for (int j =0; j < N ; j++)
    {
      if (real(fHat[j])/N > 1.e-08)
        {
          numModes +=  1;
          maxMode = j;
        }
    }
  if (numModes > 1) 
    {
      return -1;
    }
  else
    {
      return maxMode;
    }
}
int main(int argc, char* argv[])
{
  int M;
  int inputMode;
  string fft_string;
  
  cout << "input log_2(N), N = number of points" << endl;
  cin >> M ;
  cout << "input test mode < N" << endl;
  cin >> inputMode;
  cout << "input FFT being tested: Recursive, BRI, FFTW" << endl;
  cin >> fft_string;
  FFT1D* p_fft;
  FFT1DRecursive* p_fft1dR;
  FFT1DBRI* p_fft1dBRI;
  FFTW1D* p_fft1dFFTW; 

  if (fft_string == "Recursive")
    {
      p_fft1dR = new FFT1DRecursive(M);
      p_fft = dynamic_cast<FFT1D*>(p_fft1dR);
    }
  else if (fft_string == "BRI")
    {
      p_fft1dBRI = new FFT1DBRI(M);
      p_fft = dynamic_cast<FFT1D*>(p_fft1dBRI);
    }
  else if (fft_string == "FFTW")
    {
      p_fft1dFFTW = new FFTW1D(M);
      p_fft =  dynamic_cast<FFT1D*>(p_fft1dFFTW);
    }
  else
    {
      cout << "invalid FFT" << endl;
      abort();
    }
  double error = test1(p_fft);
  int mode = test2(p_fft,inputMode);
  cout << fft_string<< ": test 1: error in Gaussian  = " << error << endl;
  cout << fft_string << ": test 2: reproducing input Fourier mode " << inputMode << 
  " , output mode " << mode <<endl;
  if (fft_string == "Recursive")
    {
      delete p_fft1dR;
    }
  else if (fft_string == "BRI")
    {
      delete p_fft1dBRI;
    }
  else if (fft_string == "FFTW")
    {
      delete p_fft1dFFTW;
    }
  else
    {
      cout << "invalid FFT" << endl;
      abort();
    }
};
