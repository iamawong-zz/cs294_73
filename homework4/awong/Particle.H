#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "DIM.H"
#include "Real.H"

class DX
{
public:
  DX()
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]=0;
      }
  } 
  Real m_x[DIM];
  void increment(double a_scale, const DX& a_rhs)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]+=a_scale*a_rhs.m_x[i];
      }
  }
  void operator*=(double a_scale)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]*=a_scale;
      }
  }
};

class Particle
{
public:
  Real m_x[DIM];
  Real strength;
  void increment(const DX& a_shift)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]+=a_shift.m_x[i];
      }
  }
};

#endif
