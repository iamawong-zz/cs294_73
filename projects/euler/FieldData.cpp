#include "FieldData.H"
#include "DeltaVelocity.H"

FieldData::FieldData(){
}

//FIELD DATA is Q from the notes
FieldData::FieldData(std::tr1::shared_ptr<FFT1D> a_fft1dptr, int a_nghosts){
  m_ghosts = a_nghosts;
  m_fft1dptr = a_fft1dptr;
  m_M = a_fft1dptr->getM();
  m_N = a_fft1dptr->getN();
  int lowerCorner[2] = {0,0};
  int upperCorner[2] = {m_N-1, m_N-1};
  m_grid = Box(lowerCorner, upperCorner);
  for (int i = 0; i<DIM; i++)
    {
      m_data[i].define(m_grid.grow(a_nghosts));
    }  
//make a box with low corner = 00, high corner = m_n, m_M?? make this box m_grid
  //  m_grid = m_data.getBox();
}

const MDArray<double>& FieldData::operator[](int a_component) const{
  return m_data[a_component];  
}

MDArray<double>& FieldData::operator[](int a_component){
  return m_data[a_component];  
}

void FieldData::copyTo(FieldData& a_FieldData) const
{
  for (unsigned int i =0; i < DIM; i++)
    {
      a_FieldData[i].define(m_grid.grow(m_ghosts));
      m_data[i].copyTo(a_FieldData[i]);
    }
  a_FieldData.m_M = m_M;
  a_FieldData.m_N = m_N;
  a_FieldData.m_ghosts = m_ghosts;
  a_FieldData.m_grid = Box(m_grid);
  a_FieldData.m_fft1dptr = m_fft1dptr;
}

void FieldData::fillGhosts(MDArray<double>& array){
  int highCorner[DIM];
  int lowCorner[DIM];
  array.getBox().getHighCorner(highCorner);
  array.getBox().getLowCorner(lowCorner);
  int ghostx = lowCorner[0];
  int ghosty = lowCorner[1];
  //fill Ghosts for bottom row, top row
  for(ghostx=lowCorner[0]+1; ghostx <highCorner[0]; ghostx++){
    int bottomIndex[DIM] = {ghostx, highCorner[1]-1};
    int bottomGhostIndex[DIM] = {ghostx, lowCorner[1]};
    int topIndex[DIM] = {ghostx, lowCorner[1] + 1};
    int topGhostIndex[DIM] = {ghostx, highCorner[1]};
    array[bottomGhostIndex] = array[bottomIndex];
    array[topGhostIndex] = array[topIndex];
  }
  // fillGhosts for left row, right row
  for(ghosty=lowCorner[1]+1; ghosty < highCorner[1]; ghosty++){
    int leftIndex[DIM] = {highCorner[0] -1, ghosty};
    int leftGhostIndex[DIM] = {lowCorner[0], ghosty};
    int rightIndex[DIM] = {lowCorner[0]+1, ghosty};
    int rightGhostIndex[DIM] = {highCorner[0], ghosty};
    array[leftGhostIndex] = array[leftIndex];
    array[rightGhostIndex] = array[rightIndex];
    }
}

void FieldData::fillGhosts(){
  fillGhosts(m_data[0]);
  fillGhosts(m_data[1]);
  // DeltaVelocities do not have the ghost data
  // Fill in the edges from the mirror on the opposite side of the box
  //Does this mean that I want to expand m_grid to have extra spaces around the outside?
}


void FieldData::increment(const double& a_scalar,
			  const DeltaVelocity& a_fieldIncrement){
  Box bi =a_fieldIncrement.m_data[0].getBox();
  MDArray<Real> temp(bi);
  for (unsigned int i = 0; i <DIM; i++){
    a_fieldIncrement.m_data[i].copyTo(temp);
    temp*=a_scalar;
    m_data[i] += temp;
  }
}

// (ParticleShift& a_kOut, 
//                      const Real& a_time, const Real& dt, 
//                      const ParticleSet& a_state,
//                      const ParticleShift& a_kIn)

