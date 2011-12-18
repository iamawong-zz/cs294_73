# include "Box.H"
#include <cassert>
Box::Box()
{
  for (int dir = 0 ; dir < DIM ; dir++)
    {
      m_lowCorner[dir] = 0;     
      m_highCorner[dir] = -1;
    }
  m_isEmpty = true;
}
Box::Box(const Box& a_box)
{
  *this = a_box;
}
Box::Box(const int a_lowCorner[DIM],const int a_highCorner[DIM])
{
  m_isEmpty = false;
  for (int i = 0 ; i < DIM; i++)
    {
      m_isEmpty = m_isEmpty || (a_lowCorner[i] > a_highCorner[i]);
    }
  if (m_isEmpty)
    for (int dir = 0 ; dir < DIM ; dir++)
    {
      m_lowCorner[dir] = 0;     
      m_highCorner[dir] = -1;
    }
  else
    {   
      for (int dir = 0 ; dir < DIM ; dir++)
        {
          m_lowCorner[dir] = a_lowCorner[dir];
          m_highCorner[dir] = a_highCorner[dir];
        }
    }
}
Box Box::operator*(const Box& a_rightBox) const
{
  int newLo[DIM];
  int newHi[DIM];
  if (isEmpty() || a_rightBox.isEmpty())
    {
      return Box();
    }
  bool emptyFlag = false;
  for (int i = 0 ; i < DIM; i++)
    {
      newLo[i] = m_lowCorner[i];
      if (m_lowCorner[i] < a_rightBox.m_lowCorner[i])
        {
          newLo[i] = a_rightBox.m_lowCorner[i];
        }
      newHi[i] = m_highCorner[i];
      if (m_highCorner[i] > a_rightBox.m_highCorner[i])
        {
          newHi[i] = a_rightBox.m_highCorner[i];
        }
      emptyFlag = emptyFlag || (newHi[i] < newLo[i]);
    }
  if (emptyFlag)
    {
      return Box();
    }
  else
    {
      return Box(newLo,newHi);
    }
}
Box Box::shift(int a_tuple[DIM]) const
{
  Box returnBox = Box(*this);
  if (m_isEmpty) {return returnBox;}
  for (int dir = 0 ; dir < DIM ; dir++)
    {
      returnBox.m_lowCorner[dir] += a_tuple[dir];
      returnBox.m_highCorner[dir] += a_tuple[dir];
    }
  return returnBox;
}
Box Box::shift(int a_direction, int a_offset) const
{
  Box returnBox = Box(*this);
  if (m_isEmpty) {return returnBox;}
  returnBox.m_lowCorner[a_direction] += a_offset;
  returnBox.m_highCorner[a_direction] += a_offset;
  return returnBox;
}
Box Box::grow(int a_offset) const
{
  Box returnBox(*this);
  for (int i=0; i< DIM;i++)
    {
      returnBox.m_lowCorner[i] -= a_offset;
      returnBox.m_highCorner[i] += a_offset;
    }
  return returnBox;
}
Box Box::grow(int a_dir, int a_offset) const
{
  Box returnBox(*this);
 
  returnBox.m_lowCorner[a_dir] -= a_offset;
  returnBox.m_highCorner[a_dir] += a_offset;
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
  if (m_isEmpty)
    {
      return 0;
    }
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
int Box::length(int a_dir) const
{
  return m_highCorner[a_dir] - m_lowCorner[a_dir] + 1;
}
bool Box::isSquare() const
{
  int l0 = length(0);
  bool flag = true;
  for (int dir = 0; dir < DIM; dir++)
    {
      flag = flag&&(l0 == length(dir));
    }
  return flag;
}
const bool& Box::isEmpty() const
{
  return m_isEmpty;
}

