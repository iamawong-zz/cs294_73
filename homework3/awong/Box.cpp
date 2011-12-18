# include "Box.H"
#include <cassert>
Box::Box()
{
  for (int dir = 0 ; dir < DIM ; dir++)
    {
      m_lowCorner[dir] = 0;     
      m_highCorner[dir] = -1;
    }
}
Box::Box(const Box& a_box)
{
  for (int dir = 0 ; dir < DIM ; dir++)
    {
      m_lowCorner[dir] = a_box.m_lowCorner[dir];
      m_highCorner[dir] = a_box.m_highCorner[dir];
    }
}
Box::Box(const int a_lowCorner[DIM],const int a_highCorner[DIM])
{
  for (int dir = 0 ; dir < DIM ; dir++)
    {
      m_lowCorner[dir] = a_lowCorner[dir];
      m_highCorner[dir] = a_highCorner[dir];
    }
}
//Box returns a new box, should have put const to make clear.
Box Box::operator*(const Box& a_rightBox) const
{
  int newLo[DIM];
  int newHi[DIM];
  for (int i = 0 ; i < DIM; i++)
    {
      newLo[i] = m_lowCorner[i];
      if (m_lowCorner[i] < a_rightBox.m_lowCorner[i])
        {
          newLo[i] = a_rightBox.m_lowCorner[i];
        }
      newLo[i] = m_highCorner[i];
      if (m_highCorner[i] > a_rightBox.m_highCorner[i])
        {
          newHi[i] = a_rightBox.m_highCorner[i];
        }
    }
  return Box(newLo,newHi);
}
Box Box::shift(int a_direction, int a_offset)
{
  Box returnBox = Box(*this);
  returnBox.m_lowCorner[a_direction] += a_offset;
  returnBox.m_highCorner[a_direction] += a_offset;
  return returnBox;
}
Box Box::grow(int a_offset)
{
  Box returnBox(*this);
  for (int i=0; i< DIM;i++)
    {
      returnBox.m_lowCorner[i] -= a_offset;
      returnBox.m_highCorner[i] += a_offset;
    }
  return returnBox;
}
void Box::getLowCorner(int a_lowcorner[DIM]) const
{
  for (int i = 0;i < DIM; i++)
    {
      a_lowcorner[i] = m_lowCorner[i];
    }
}
void Box::getHighCorner(int a_highcorner[DIM]) const
{
  for (int i = 0;i < DIM; i++)
    {
       a_highcorner[i] = m_highCorner[i];
    }
}
int Box::linearIndex(const int a_tupleIndex[DIM]) const
{
  int factor = 1;
  int linIndex = 0;
  for (int i = 0;i < DIM;i++)
    { 
      linIndex = linIndex + (a_tupleIndex[i] - m_lowCorner[i])*factor;
      factor = factor*(m_highCorner[i] - m_lowCorner[i]+1);
    }
  return linIndex;
}
     
void Box::tupleIndex(int a_linearIndex, int a_tupleIndex[DIM]) const
{
  int size = sizeOf();
  for (int i=DIM-1 ; i >=1 ; i--)
    {
      size = size/(m_highCorner[i] - m_lowCorner[i] + 1);
      int remainder = a_linearIndex%size;
      a_tupleIndex[i] =  (a_linearIndex-remainder)/size + m_lowCorner[i];
      a_linearIndex = remainder;
      
    }
  a_tupleIndex[0] = a_linearIndex + m_lowCorner[0];
}
int Box::sizeOf() const
{
  int retval = 1;
  for (int dir = 0;dir < DIM;dir++)
    {
      retval=(m_highCorner[dir] - m_lowCorner[dir] + 1)*retval;
    }
  return retval;
}
bool Box::operator==(const Box& a_rhsBox) const
{
  bool retval = true;
  for (int dir=0 ; dir < DIM ; dir++)
    {
      retval = retval && (m_lowCorner[dir] == a_rhsBox.m_lowCorner[dir]);
      retval = retval && (m_highCorner[dir] == a_rhsBox.m_highCorner[dir]);
    }
  return retval;
}
