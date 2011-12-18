#include "ParticleSet.H"
#include <cassert>
#include <iostream>

void ParticleShift::increment(double a_scale, const ParticleShift& a_rhs)
{
  assert(m_particles.size() == a_rhs.m_particles.size());
  for(unsigned int i=0; i<m_particles.size(); i++)
    {
      m_particles[i].increment(a_scale, a_rhs.m_particles[i]);
    }
}

void ParticleShift::operator*=(double a_scale)
{
  for(unsigned int i=0; i<m_particles.size(); i++)
    {
      m_particles[i]*=a_scale;
    }
}
  
void ParticleShift::setToZero()
{
  for(unsigned int i=0; i<m_particles.size(); i++)
    {
      m_particles[i]*=0;
    }
}

void ParticleShift::init(const ParticleSet& a_rhs)
{
  m_particles.resize(0);
  m_particles.resize(a_rhs.m_particles.size());
}

void ParticleSet::increment(const ParticleShift& a_rhs)
{
  assert(m_particles.size() == a_rhs.m_particles.size());
  for(unsigned int i=0; i<m_particles.size(); i++)
    {
      m_particles[i].increment(a_rhs.m_particles[i]);
    }
}
void ParticleSet::rebin(const Box& a_domain, Real a_dx, const Real a_lowCorner[DIM])
{
  int locationTuple[DIM];
  m_dx=a_dx;
  m_box=a_domain;
  (*m_bins).define(m_box);
  
  for(unsigned int k=0; k<m_particles.size(); k++)
    {
      locationTuple[0] = floor((m_particles[k].m_x[0]-a_lowCorner[0])/a_dx + 0.5);
      locationTuple[1] = floor((m_particles[k].m_x[1]-a_lowCorner[0])/a_dx + 0.5);
      (*m_bins)[locationTuple].push_back(k);
    }
}
