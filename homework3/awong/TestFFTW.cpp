
#include <fftw3.h>
#include <cmath>

void fillInSignal(fftw_complex* a_in, int a_n, int terms, double* a_frequencies);

int main(int argc, char** argv)
{
  int N = 2048;  double f[4] = {1.02, 2.0, 4.002, 8.01};
  fftw_complex *in, *out;
  fftw_plan p;
        
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  for(int i=0; i<4; i++)
    {
      fillInSignal(in, N, i, f);
      fftw_execute(p);
    }
  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}


void fillInSignal(fftw_complex* a_in, int a_n, int terms, double* a_frequency)
{
  // cheat the system and put in a real sequence of values into
  // the complex number as the even and odd points
  double x=0;
  double dx = 2*M_PI/(a_n-1);
  for(int i=0; i<a_n; i++)
    {
      a_in[0][0] = sin(x);
      a_in[0][1] = sin(x+dx/2.0);
      a_in++;
      x+=dx;
    }
}
