#include "FFTW1DSine.H"
#include "fftw3.h"
#include <cstdlib>

using namespace std;
FFTW1DSine::FFTW1DSine(){};
FFTW1DSine::FFTW1DSine(unsigned int a_N)
{
  m_N = a_N;
  m_in.resize(m_N);
  m_out.resize(m_N);

  m_plan= fftw_plan_r2r_1d( m_N, &m_in[0], &m_out[0],
                                      FFTW_RODFT00, FFTW_ESTIMATE);
}

FFTW1DSine::~FFTW1DSine()
{
  fftw_destroy_plan(m_plan);
}

void FFTW1DSine::transform(vector<double > & a_fHat, 
                           const vector<double> & a_f) const
{
  for (unsigned int k = 0;k < m_in.size();k++)
    {
      m_in[k] = a_f[k];
    }
  fftw_execute(m_plan);
 for (unsigned int k = 0;k < m_in.size();k++)
    {
      a_fHat[k] = m_out[k];
    }
}
