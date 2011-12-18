#include "FFTW1D.H"

FFTW1D::FFTW1D(unsigned int a_M)
{
  /// Sizing everything correctly to 2^m_M
  m_M = a_M;
  m_N = Power(2,m_M);
  m_in.resize(m_N);
  m_out.resize(m_N);

  /// Make the plans in the constructor
  fftw_complex *in;
  fftw_complex *out;
  in = reinterpret_cast<fftw_complex*>(&(m_in[0]));
  out = reinterpret_cast<fftw_complex*>(&(m_out[0]));

  m_forward = fftw_plan_dft_1d(m_N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  m_inverse = fftw_plan_dft_1d(m_N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
}

void FFTW1D::forwardFFTCC(vector<complex<double> > & a_fHat,
			  const vector<complex<double> >& f) const
{
  /*  fftw_complex *in;
  fftw_complex *out;
  unsigned int N = f.size();
  fftw_plan p;

  /// redirecting it because it's const
  m_in = f;

  /// recast the input to fftw_complex
  in = reinterpret_cast<fftw_complex*>(&(m_in[0]));
  out = reinterpret_cast<fftw_complex*>(&(a_fHat[0]));

  p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);*/

  m_in = f;

  fftw_execute(m_forward); /* repeat as needed */

  a_fHat = m_out;

  //fftw_destroy_plan(m_forward);
  //fftw_free(in); fftw_free(out);  
}

void FFTW1D::inverseFFTCC(vector<complex<double> > & a_f,                                                                                            
                          const vector<complex<double> > & a_fHat) const
{
  /*fftw_complex *in;
  fftw_complex *out;
  fftw_plan p;
  unsigned int N = a_fHat.size();

  m_in = a_fHat;

  in = reinterpret_cast<fftw_complex*>(&(m_in[0]));
  out = reinterpret_cast<fftw_complex*>(&(a_f[0]));

  p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  */
  
  m_in = a_fHat;
  
  fftw_execute(m_inverse); /* repeat as needed */
  
  a_f = m_out;

  //fftw_destroy_plan(m_inverse);
  //fftw_free(in); fftw_free(out);
}

/// Destructor and destroy the plans
FFTW1D::~FFTW1D()
{
  fftw_destroy_plan(m_forward);
  fftw_destroy_plan(m_inverse);
}
