#ifndef _PARTICLEVELOCITIES_H_
#define _PARTICLEVELOCITIES_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "MDArray.H"
#include "FFTW1DSine.H"
#include "FFTMDSine.H"
#include "PoissonInfinite.H"
#include "ParticleSet.H"
#include "MLCVortexMethod.H"
using namespace std;
class ParticleVelocities
{
public:
  void operator()(ParticleShift& a_kOut, 
                     double a_time, double dt, 
                     const ParticleSet& a_state,
                     const ParticleShift& a_kIn);
};
#endif
