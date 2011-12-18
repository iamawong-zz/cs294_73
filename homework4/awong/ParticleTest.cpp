

#include "RK4.H"
#include <vector>
#include <cmath>
#include <iostream>
#include "ParticleSet.H"
#include "WriteMDArray.H"
#include <cmath>

Real omega = 0.10;
class F
{
public:
  void operator()(ParticleShift& a_result, double a_time, double a_dt, const ParticleSet& a_val, const ParticleShift& a_shift)
  {
    double rotationOrigin[DIM] = {0.5, 0.5};
    double qCurrentState[DIM];
    assert(a_result.m_particles.size() == a_val.m_particles.size());
    assert(a_result.m_particles.size() == a_shift.m_particles.size());

    for (unsigned int i=0; i<a_result.m_particles.size(); i++)
      {
	qCurrentState[0] = a_val.m_particles[i].m_x[0]+a_shift.m_particles[i].m_x[0];
	qCurrentState[1] = a_val.m_particles[i].m_x[1]+a_shift.m_particles[i].m_x[1];
	a_result.m_particles[i].m_x[0] = -a_dt*(qCurrentState[1]-rotationOrigin[1])*omega;
	a_result.m_particles[i].m_x[1] = a_dt*(qCurrentState[0]-rotationOrigin[0])*omega;
      }
  }
  
 
};


int main(int argc, char* argv[])
{
  
  RK4<ParticleSet, F, ParticleShift> integrator;
  
  ParticleSet initialCondition;
  initialCondition.m_particles.resize(100);
  for(int i=0; i<20; i++)
    {
      Particle& p = initialCondition.m_particles[i];
      for(int d=0; d<DIM; d++)
	{
	  p.m_x[d]=drand48();
       
	}
    }

 

  std::vector<Real> error;
  std::vector<Real> timestep;
  std::vector<ParticleSet> evolutions;
  evolutions.reserve(10);
  for(int m=10; m<=160; m*=2)
    {

      double time = 0;
      evolutions.push_back(initialCondition);

      double dt = 2*M_PI/(m*omega);
      timestep.push_back(dt);
      char name[20];
      int i=0;
      for(; i<m; i++)
	{
	  sprintf(name, "Part.%d.%d",m,i);
	  PWrite(name, &(evolutions.back()));
	  integrator.advance(time, dt, evolutions.back());
	  time+=dt;
	}
      sprintf(name, "Part.%d.%d",m,i);
      PWrite(name, &(evolutions.back()));

    }
  for(unsigned int i=0; i<evolutions.size()-1; i++)
    {
      ParticleSet& coarse = evolutions[i];
      ParticleSet& fine = initialCondition;
      double sumS=0;
      for(unsigned int p=0; p<coarse.m_particles.size(); p++)
	{
	  Particle& c = coarse.m_particles[p];
	  Particle& f = fine.m_particles[p];
	  for(int d=0; d<DIM; d++)
	    {
	      double m = c.m_x[d]-f.m_x[d];
	      sumS += m*m;
	    }
	}
      //sumS*=timestep[i];
      double e = sqrt(sumS);
      std::cout<<"L2 norm error in evolution dt="<<timestep[i]<<" : "<<e<<std::endl;
      error.push_back(e);

    }
  for(unsigned int i=0; i<error.size()-1; i++)
    {
      double base = timestep[i]/timestep[i+1];
      double errorOrder = log(error[i]/error[i+1])/log(base);
      std::cout<<"order = "<<errorOrder;
      std::cout<<std::endl;
    }
  
  
}
