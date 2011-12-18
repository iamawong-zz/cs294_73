#ifndef _PARTICLESET_H_
#define _PARTICLESET_H_

#include <vector>
#include "Particle.H"
#include "MDArray.H"
#include <list>
#include <tr1/memory>

using namespace std;
typedef MDArray<list<unsigned int> > Bins;

class ParticleSet;

class ParticleShift
{
public:
  vector<DX> m_particles;
  void init(const ParticleSet& a_particles);
  void increment(double a_scale, const ParticleShift& a_rhs);
  void operator*=(double a_scale);
  void setToZero();
};


class ParticleSet
{
public:
  vector<Particle> m_particles;
  std::tr1::shared_ptr<Bins> m_bins;
  Real m_dx;
  Real m_cutoff;
  Box m_box;
  int m_corrRadius;
  Real m_lowCorner[DIM];

  void increment(const ParticleShift& a_shift);
  void rebin(const Box& a_domain, Real a_dx, const Real a_lowCorner[DIM]);
};


#endif
