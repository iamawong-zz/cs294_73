#include "DeltaVelocity.H"
#include "FieldData.H"


DeltaVelocity::DeltaVelocity()
{
  // do nothing
}

/* DeltaVelocity::~DeltaVelocity()
{
  // do nothing
  }*/


void DeltaVelocity::init(const FieldData& a_vel){
  m_grid = a_vel.m_grid.grow(-1); 
 for(int i =0; i < DIM; i++)
    {
      // the other stuff doesn't work. fixing it ;-)
      // -Johnny
      m_data[i].define(a_vel.m_data[i].getBox().grow(-1));
      m_data[i].setVal(0.);
    }
  m_grid = a_vel.m_grid;
}


void DeltaVelocity::operator*=(double a_scalar)
{
  for(unsigned int i = 0; i < DIM; i++)
    {
      m_data[i]*=a_scalar;
    }
}

const MDArray<double>& DeltaVelocity::operator[](int a_component) const{
  return m_data[a_component];
}

MDArray<double>& DeltaVelocity::operator[](int a_component) {
  return m_data[a_component];
}


void DeltaVelocity::increment(const double& a_scalar,
			      const DeltaVelocity& a_fieldIncrement){
  Box bi =a_fieldIncrement.m_data[0].getBox();
  MDArray<Real> temp(bi);
  for (unsigned int i = 0; i <DIM; i++){
    a_fieldIncrement.m_data[i].copyTo(temp);
    temp*=a_scalar;
    m_data[i] += temp;
  }
 }

//use the MDArray += and *= in order to do this, it will do all of this on the valid region
