#include "FFT1DRecursive.H"

FFT1DRecursive::FFT1DRecursive(const unsigned int& a_M):FFT1D(a_M)
{
};

/// Forward FFT: a_fHat[k] = \sum_j=0^{N-1} a_f[j] z^{j k}, z = e^{-2 \pi \iota /m_N}                                                                        
void FFT1DRecursive::forwardFFTCC(vector<complex<double> > & a_fHat,
			  const vector<complex<double> >& f) const
{  
  unsigned int currLength = f.size();
  /// Compute the N=2 FFT
  if(currLength == 2)
    {
      a_fHat[0] = f[0] + f[1];
      a_fHat[1] = f[0] - f[1];
    }
  else
    {
      /// Initialize some vectors to hold stuff
      vector<complex<double> > f_even(currLength/2, 0);
      vector<complex<double> > f_odd(currLength/2, 0);
      vector<complex<double> > fHat_even(currLength/2, 0);
      vector<complex<double> > fHat_odd(currLength/2, 0);
      
      /// Loops to split up the even and odd indicies
      for(unsigned int i=0; i<currLength/2; i++)
	{
	  f_even[i] = f[2*i];
	  f_odd[i] = f[2*i+1];
	}
  
      /// RECURSION is fun. But just do FFT on these bad boys.
      forwardFFTCC(fHat_even, f_even);
      forwardFFTCC(fHat_odd, f_odd);

      /// Piece back together after FFT on even and odd indicies
      for(unsigned int k=0; k<currLength/2; k++)
	{
	  double thetaArg = -2*M_PI*k/currLength;
	  complex<double> z (cos(thetaArg), sin(thetaArg));
	  a_fHat[k] = fHat_even[k] + z*fHat_odd[k];
	  a_fHat[k+currLength/2] = fHat_even[k] - z*fHat_odd[k];
	}
    }
  return;
};
  
/// inverse FFT: a_f[j] = \sum_{k=0}^{N-1} a_fHat[k] z^{j k}, z = e^{2 \pi \iota /m_N}                                                                       
void FFT1DRecursive::inverseFFTCC(vector<complex<double> > & a_f,
			  const vector<complex<double> > & a_fHat) const
{  
  unsigned int currLength = a_fHat.size();
  /// Compute the N=2 FFT
  if(currLength == 2)
    {
      a_f[0] = a_fHat[0] + a_fHat[1];
      a_f[1] = a_fHat[0] - a_fHat[1];
    }
  else
    {
      /// Initialize some vectors to hold stuff
      vector<complex<double> > fHat_even(currLength/2, 0);
      vector<complex<double> > fHat_odd(currLength/2, 0);
      vector<complex<double> > f_even(currLength/2, 0);
      vector<complex<double> > f_odd(currLength/2, 0);

      /// Loops to split up the even and odd indicies
      for(unsigned int i=0; i<currLength/2; i++)
	{
	  fHat_even[i] = a_fHat[2*i];
	  fHat_odd[i] = a_fHat[2*i+1];
	}
  
      /// RECURSION is fun. But just do FFT on these bad boys.
      inverseFFTCC(f_even, fHat_even);
      inverseFFTCC(f_odd, fHat_odd);

      /// Piece back together after FFT on even and odd indicies
      for(unsigned int k=0; k<currLength/2; k++)
	{
	  double thetaArg = 2*M_PI*k/currLength;
	  complex<double> z (cos(thetaArg), sin(thetaArg));
	  a_f[k] = f_even[k] + z*f_odd[k];
	  a_f[k+currLength/2] = f_even[k] - z*f_odd[k];
	}
    }
  return;
};
