#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <complex>
#include "Box.H"
#include "MDArray.H"
#include "Real.H"
#include "MLCVortexMethod.H"
#include "WriteMDArray.H" 
#include "VisitWriter.H"
using namespace std;
MLCVortexMethod::MLCVortexMethod(const Real& a_cutoff, 
                                 const Real& a_h, const Real a_lowCorner[DIM])
{
  m_cutoff = a_cutoff;
  m_h = a_h;
  m_lowCorner[0] = a_lowCorner[0];
  m_lowCorner[1] = a_lowCorner[1];

  int stencilOffset[MLCSTENCILSIZE][DIM] = {{0,0},{1,0},{-1,0},{0,1},{0,-1}};
    
  for (int l = 0; l < MLCSTENCILSIZE;l++)
    {
      for (int dir = 0;dir < DIM;dir++)
        {
          m_stencilOffset[l][dir] = stencilOffset[l][dir];
        }
      m_zpoints[l] = complex<double>(stencilOffset[l][0],
                                     stencilOffset[l][1]);
      // cout << m_zpoints[l] << endl;
    }
}
void 
MLCVortexMethod::computeVelocities(MDArray<Real> a_gridVelocities[DIM],
                                   const ParticleSet& a_p,
                                   int a_tuple[DIM]) const
{
  int tuple[DIM];
  Real vel[DIM];
  int tupleDest[DIM];
  
  a_gridVelocities[0].setVal(0.);
  a_gridVelocities[1].setVal(0.);
  
  for (int kDest = 0; kDest < a_gridVelocities[0].getBox().sizeOf();kDest++)
    {
      a_gridVelocities[0].getBox().tupleIndex(kDest,tupleDest);
      //for (int k = 0;k < a_subBox.sizeOf();k++)
        {
          // a_subBox.tupleIndex(k,tuple);
          
          list<unsigned int>::iterator it;
          std::tr1::shared_ptr<Bins> bins=a_p.m_bins;
          for (it = (*bins)[a_tuple].begin();it != (*bins)[a_tuple].end();it++)
            {
              getVel(vel,tupleDest,a_p.m_particles[*it]);
              a_gridVelocities[0][kDest] += vel[0];
              a_gridVelocities[1][kDest] += vel[1];
            }
        }
    }
};
void 
MLCVortexMethod::MLCInterpolate(ParticleShift& a_kOut,
                                const ParticleSet& a_p,
                                const MDArray<Real> a_gridVelocities[DIM],
                                int a_bin[DIM]) const
{
  int tuple[DIM];
  Real velStencil[MLCSTENCILSIZE][DIM];
  Box bx(a_bin,a_bin);
  Box localCells = bx.grow(a_p.m_corrRadius);

  for (int l = 0; l < MLCSTENCILSIZE;l++)
    {
      
      for (int dir =0;dir < DIM; dir++)
        {
          tuple[dir] = a_bin[dir] + m_stencilOffset[l][dir];
        }
      velStencil[l][0] = a_gridVelocities[0][tuple];
      velStencil[l][1] = a_gridVelocities[1][tuple];
    }

  for (int k=0; k < localCells.sizeOf();k++)
    {
      localCells.tupleIndex(k,tuple);
      list<unsigned int>::iterator it;
      std::tr1::shared_ptr<Bins> bins=a_p.m_bins;
      for (it = (*bins)[tuple].begin();it != (*bins)[tuple].end();it++)
        {
          int sourceIndex = *it;
          
          for (int l = 0; l < MLCSTENCILSIZE;l++)
            {
              int tupleDest[DIM];
              for (int dir = 0;dir < DIM;dir++)
                {
                  tupleDest[dir] = a_bin[dir] + m_stencilOffset[l][dir];
                }
              Real vel[DIM];
              getVel(vel,tupleDest,a_p.m_particles[sourceIndex]);
               velStencil[l][0] -= vel[0];
               velStencil[l][1] -= vel[1];
            }
          list<unsigned int>::iterator it2;
          for (it2 = (*bins)[a_bin].begin();it2 != (*bins)[a_bin].end();it2++)
            {
              int destIndex = *it2;
              Real vel[DIM];
              
              getVel(vel,a_p.m_particles[destIndex],
                     a_p.m_particles[sourceIndex]);
               a_kOut.m_particles[destIndex].m_x[0] += vel[0];
               a_kOut.m_particles[destIndex].m_x[1] += vel[1];
            }
        }      
    }
  interpolate(a_kOut,a_p,velStencil,a_bin);
};
void MLCVortexMethod::interpolate(ParticleShift& a_kOut,
                                  const ParticleSet& a_p,
                                  const Real a_velStencil[MLCSTENCILSIZE][DIM],
                                  int a_tuple[DIM]) const
{
  // create complex stencil.
  complex<Real> a[MLCSTENCILSIZE],coef[MLCSTENCILSIZE];
  complex<Real> z0 = complex<Real> (a_tuple[0],a_tuple[1]);
  complex<double> ii(0,1);
  complex<double> four(4,0);
  for (int l = 0;l < MLCSTENCILSIZE;l++)
    {
      a[l] = complex<Real> (a_velStencil[l][0],-a_velStencil[l][1])/four;
    }
  
  
  coef[0] = four*a[0];
  coef[1] = a[1] - a[2] - ii*a[3] + ii*a[4];
  coef[2] = a[1] + a[2] - a[3] - a[4];
  coef[3] = a[1] - a[2] + ii*a[3] - ii*a[4];
  coef[4] = -four*a[0] + a[1] + a[2] + a[3] +a[4];
  
  list<unsigned int>::iterator it;
  std::tr1::shared_ptr<Bins> bins=a_p.m_bins;
  for (it = (*bins)[a_tuple].begin();it != (*bins)[a_tuple].end();it++)
    {
      const Particle& particle = a_p.m_particles[*it];
      complex<Real> zAbs = complex<Real>(particle.m_x[0]/m_h,particle.m_x[1]/m_h);
      complex<Real> z = (zAbs - z0);
      complex<Real> cvelocity = 
        (((coef[4]*z + coef[3])*z + coef[2])*z + coef[1])*z + coef[0];
      // cout << zAbs << " " <<  cvelocity << endl;
      // cout << a_kOut.m_particles[*it].m_x[0] << endl;
      // cout << a_kOut.m_particles[*it].m_x[1] << endl;
      a_kOut.m_particles[*it].m_x[0] += real(cvelocity);
      a_kOut.m_particles[*it].m_x[1] -= imag(cvelocity);
    }
};
void MLCVortexMethod::getVel(Real a_vel[DIM], 
                             int a_tuple[DIM],
                             const Particle& a_source) const
{
  Real dx = (m_lowCorner[0] + a_tuple[0])*m_h - a_source.m_x[0];
  Real dy = (m_lowCorner[1] + a_tuple[1])*m_h - a_source.m_x[1];
  Real str = a_source.strength;
  Real r = sqrt(pow(dx,2) + pow(dy,2));
  Real rcut = max(r,m_cutoff);
  if (r > 1.e-12*m_h)
    {
      a_vel[0] = -str*dy/r/rcut/(2*M_PI);
      a_vel[1] =  str*dx/r/rcut/(2*M_PI);
    }
    else
    {
      a_vel[0] = 0.;
      a_vel[1] = 0.;
    }
};
void MLCVortexMethod::getVel(Real a_vel[DIM], 
                             const Particle& a_dest,
                             const Particle& a_source) const

{
  Real dx = a_dest.m_x[0] - a_source.m_x[0];
  Real dy = a_dest.m_x[1] - a_source.m_x[1];
  Real str = a_source.strength;
  Real r = sqrt(pow(dx,2) + pow(dy,2));
  Real rcut = max(r,m_cutoff);
   if (r > 1.e-12*m_h)
    {
      a_vel[0] = -str*dy/r/rcut/(2*M_PI);
      a_vel[1] =  str*dx/r/rcut/(2*M_PI);
    }
    else
    {
      a_vel[0] = 0.;
      a_vel[1] = 0.;
    }
};

