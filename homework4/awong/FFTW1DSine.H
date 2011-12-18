#ifndef _FFTW1DSINE_H_
#define _FFTW1DSINE_H_
#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include "PowerItoI.H"
#include "fftw3.h"

using namespace std;
class FFTW1DSine
{
public:
  FFTW1DSine();
  ~FFTW1DSine();
  FFTW1DSine(unsigned int a_N);
  void transform(vector<double>  & a_fHat, 
                         const vector<double> & f) const;
  

protected:
  mutable vector<double>  m_in;
  vector<double> m_out;
  fftw_plan m_plan;
  int m_N;
};
#endif
