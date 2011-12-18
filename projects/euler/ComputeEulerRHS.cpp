#include "ComputeEulerRHS.H"
#include "FFTW1D.H"
#include "WriteMDArray.H"

void ComputeEulerRHS::operator() (DeltaVelocity& a_newDv,
				  const double& a_time, 
				  const double& a_dt,
				  const FieldData& a_velocity,
				  DeltaVelocity& a_oldDv) {
  double a_scalar = 1.0;
  
  //std::tr1::shared_ptr<FFT1D> fft1dptr(a_velocity.m_N);  std::tr1::shared_ptr<FFT1D> fft1dptr(new FFTW1D(a_velocity.m_M));

  //create dummy variables
  FieldData velTemp;
  a_velocity.copyTo(velTemp);
  DeltaVelocity deltaV;
  //int index[DIM];
  
  Projection projVel(a_velocity.m_fft1dptr);
  Box bx = a_velocity.m_grid;
  
  //increment veltemp by a_oldDv
  velTemp.increment(a_scalar,a_oldDv);
  
  //apply advection on veltemp (what is uu?) advectionOperator(DeltaVelocity& a_divuu,const FieldData& a_velocity,const Box m_grid, const double& a_h);
  advectionOperator(deltaV,velTemp, a_velocity.m_grid,a_dt);
  velTemp.increment(-a_dt,deltaV);
  
  projVel.applyProjection(velTemp);

  //assigning a_newDv :: veltemp2 - a_velocity
  for(int i=0; i<bx.sizeOf(); i++) {
    int tuple[DIM];
    bx.tupleIndex(i, tuple);
    /*  for(int i = 0; i <a_velocity.m_M; i++) {
	for(int j = 0; j <a_velocity.m_N; j++) {
	index[0]= i;
	index[1] = j;
	a_newDv[0][index]= (velTemp2[0])[index] - (a_velocity[0])[index];
	a_newDv[1][index]= (velTemp2[1])[index] - (a_velocity[1])[index];
	}
	}*/
    a_newDv[0][tuple]= (velTemp[0])[tuple] - (a_velocity[0])[tuple];
    a_newDv[1][tuple]= (velTemp[1])[tuple] - (a_velocity[1])[tuple];
  }
  // MDWrite(a_newDv[1]);
}
