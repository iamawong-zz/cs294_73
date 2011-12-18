#ifndef _MLC_H_
#define _MLC_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <complex>
#include "Box.H"
#include "MDArray.H"
#include "ParticleSet.H"
#include "Real.H"
#define MLCSTENCILSIZE 5
class MLCVortexMethod
{
public:
  /// Constructor.
  MLCVortexMethod(const Real& a_cutoff, 
                  const Real& a_h,
                  const Real a_lowCorner[DIM]);
  /// compute the velocities on a grid induced by the particles in a_p contained
  /// in bin a_bin, and store in a_gridVelocities.
  void computeVelocities(MDArray<Real> a_gridVelocities[DIM],
                         const ParticleSet& a_p,
                         int a_bin[DIM]) const;
  /// compute the velocities induced on the particles inthe bin a_bin, 
  /// and store in a_velocities.
  void MLCInterpolate(ParticleShift& a_velocities,
                      const ParticleSet& a_particles,
                      const MDArray<Real> a_gridVelocities[DIM],
                      int a_bin[DIM]) const; 
  
private:
  void interpolate(ParticleShift& a_kOut,
                   const ParticleSet& a_p,
                   const Real a_velStencil[MLCSTENCILSIZE][DIM] ,
                   int a_tuple[DIM]) const;  
	
  void getVel(Real a_vel[DIM], 
              const Particle& a_source,
              const Particle& a_dest) const;
	
  void getVel(Real a_vel[DIM], 
              int a_dest[DIM],
              const Particle& a_source) const;
	
  Real m_cutoff;
  Real m_h;
  int m_corr;
  Real m_lowCorner[DIM];
  int m_stencilOffset[MLCSTENCILSIZE][DIM];
  complex<Real > m_zpoints[MLCSTENCILSIZE];
};
#endif
