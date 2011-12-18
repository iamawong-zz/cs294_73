#ifndef _FFTRECURSIVE_H_
#define _FFTRECURSIVE_H_
#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
using namespace std;
void FFTRecursive(vector<complex<double> >& fHat,
                  const vector<complex<double> >& f,
                  const vector<complex<double> >& a_zLevel,
                  const int& a_level);
#endif
