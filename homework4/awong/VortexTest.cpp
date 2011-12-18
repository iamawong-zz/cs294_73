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
#include "ParticleVelocities.H"
#include "WriteMDArray.H" 
#include "RK4.H"

void timeStepping(int steps, double dt, ParticleSet ps) {
  RK4<ParticleSet, ParticleVelocities, ParticleShift> integrator;
  std::vector<ParticleSet> evolutions;
  evolutions.push_back(ps);
  char name[20];
  double time=0;
  int i=0;
  for(; i<steps; i++) {
    sprintf(name, "Vort.%d.%d", steps, i);
    PWrite(name, &(evolutions.back()));
    integrator.advance(time, dt, evolutions.back());
    time+=dt;
  }
  sprintf(name, "Vort.%d.%d", steps, i);
  PWrite(name, &(evolutions.back()));
}

int main(int argc, char* argv[])
{
  unsigned int testCondition;
  cout << "input test condition that you want: " << endl;
  cin >> testCondition;
  if(testCondition == 1) {
    unsigned int N = 32;
    ParticleSet p;
    Real lowCorner[DIM];
    p.m_particles.resize(3);
    Particle& particle = p.m_particles[0];
    particle.m_x[0] = .5;
    particle.m_x[1] = .5;
    particle.strength = 1.;
    Particle& particle2 = p.m_particles[1];
    particle2.m_x[0] = .4375;
    particle2.m_x[1] = .5626;
    particle2.strength = 1.;
    Particle& particle3 = p.m_particles[2];
    particle3.m_x[0] = .45;
    particle3.m_x[1] = .55;
    particle3.strength = 1.;
    Real dx = 1./N;
    p.m_cutoff = 1.4/N;
    p.m_corrRadius = 3;
    int low[DIM] = {0,0};
    int high[DIM] = {N,N};
    Box bx(low,high);
    cout << "cutoff distance in grid units =  " << p.m_cutoff/dx << endl;
    cout << "correction radius in grid units =  " << p.m_corrRadius << endl;
    lowCorner[0] = 0.;
    lowCorner[1] = 0.;
    std::tr1::shared_ptr<Bins> bins(new Bins);
    p.m_bins = bins;
    p.rebin(bx,dx,lowCorner);
    double dt = 1.;
    int steps = 1;
    timeStepping(steps, dt, p);
  } else if(testCondition == 2) {
    unsigned int N = 32;
    ParticleSet p;
    Real lowCorner[DIM];
    p.m_particles.resize(2);
    Particle& particle = p.m_particles[0];
    particle.m_x[0] = .5;
    particle.m_x[1] = .5;
    particle.strength = 1.;
    Particle& particle2 = p.m_particles[1];
    particle2.m_x[0] = .5;
    particle2.m_x[1] = .25;
    particle2.strength = 0.;
    Real dx = 1./N;
    p.m_cutoff = 1.4/N;
    p.m_corrRadius = 3;
    int low[DIM] = {0,0};
    int high[DIM] = {N,N};
    Box bx(low,high);
    cout << "cutoff distance in grid units =  " << p.m_cutoff/dx << endl;
    cout << "correction radius in grid units =  " << p.m_corrRadius << endl;
    lowCorner[0] = 0.;
    lowCorner[1] = 0.;
    std::tr1::shared_ptr<Bins> bins(new Bins);
    p.m_bins = bins;
    p.rebin(bx,dx,lowCorner);
    double dt = 0.1;
    int steps = 300;
    timeStepping(steps, dt, p);
  } else if(testCondition == 3) {
    unsigned int N = 32;
    ParticleSet p;
    Real lowCorner[DIM];
    p.m_particles.resize(2);
    Particle& particle = p.m_particles[0];
    particle.m_x[0] = .5;
    particle.m_x[1] = .25;
    particle.strength = 1.;
    Particle& particle2 = p.m_particles[1];
    particle2.m_x[0] = .5;
    particle2.m_x[1] = .75;
    particle2.strength = 1.;
    Real dx = 1./N;
    p.m_cutoff = 1.4/N;
    p.m_corrRadius = 3;
    int low[DIM] = {0,0};
    int high[DIM] = {N,N};
    Box bx(low,high);
    cout << "cutoff distance in grid units =  " << p.m_cutoff/dx << endl;
    cout << "correction radius in grid units =  " << p.m_corrRadius << endl;
    lowCorner[0] = 0.;
    lowCorner[1] = 0.;
    std::tr1::shared_ptr<Bins> bins(new Bins);
    p.m_bins = bins;
    p.rebin(bx,dx,lowCorner);
    double dt = 0.1;
    int steps = 300;
    timeStepping(steps, dt, p);
  }
  

  return 0;
}
  
  /// and the rest of main is up to you.  build your RK4 object, implement ParticleVelocities,
  // perform time integration for the vortex problems assigned in the homework4.pdf file

