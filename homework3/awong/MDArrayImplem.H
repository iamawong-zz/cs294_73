#include "Box.H"
#include <cstdio>
#include <cmath>
#include <cassert>

template <class T> MDArray<T>::MDArray()
{
  m_isDefined = false;
  m_data = NULL;
}
template <class T> MDArray<T>::MDArray(Box a_box)
{
  m_box = a_box;
  int low[DIM];
  m_box.getLowCorner(low);
  int high[DIM];
  m_box.getHighCorner(high);
  m_data = new T[m_box.sizeOf()];
  m_isDefined=true;
}
template <class T> MDArray<T>::~MDArray()
{
  delete [] m_data;
  m_isDefined=false;
}
template <class T> void MDArray<T>::define(Box a_box)
{
  assert(m_isDefined == false);
  m_box = a_box;
  m_data = new T[m_box.sizeOf()];
  m_isDefined=true;
}
template <class T> void MDArray<T>::setVal(const T& a_val)
{
  for (int k=1;k < m_box.sizeOf();k++)
    {
      m_data[k] = a_val;
    }
}
template <class T> void MDArray<T>::operator-=(const MDArray& a_rhs)
{
  assert(m_box==a_rhs.m_box);
  for (int k=1;k < m_box.sizeOf();k++)
    {
      m_data[k] -= a_rhs.m_data[k];
    }
}
template <class T> void MDArray<T>::operator*=(const MDArray& a_rhs)
{
  assert(m_box==a_rhs.m_box);
  for (int k=1;k < m_box.sizeOf();k++)
    {
      m_data[k] *= a_rhs.m_data[k];
    }
}  
template <class T> void MDArray<T>::operator/=(const MDArray& a_rhs)
{
  assert(m_box==a_rhs.m_box);
  for (int k=0;k < m_box.sizeOf();k++)
    {
      m_data[k] /= a_rhs.m_data[k];
    }
  }      
template <class T> void MDArray<T>::operator+=(const T& a_rhs)
{
  for (int k=1;k < m_box.sizeOf();k++)
    {
      m_data[k] += a_rhs;
    }
}
template <class T> void MDArray<T>::operator-=(const T& a_rhs)
{
  for (int k=0;k < m_box.sizeOf();k++)
    {
      m_data[k] -= a_rhs;
    }
}
template <class T> void MDArray<T>::operator*=(const T& a_rhs)
{

  for (int k=0;k < m_box.sizeOf();k++)
    {
      m_data[k] *= a_rhs;
    }
}  
template <class T> void MDArray<T>::operator/=(const T& a_rhs)
{
  for (int k=0;k < m_box.sizeOf();k++)
    {
      m_data[k] /= a_rhs;
    }
}  
template <class T> MDArray<T>& MDArray<T>::shift(int a_dir,int a_len)
{
  m_box.shift(a_dir,a_len);
  return *this;
}
template <class T> T& MDArray<T>::indexShift(const int a_tuple[DIM], const int& a_dir,const int& a_shift) const
{
  int shiftedTuple[DIM];
  for (int dir = 0; dir < DIM; dir++)
    {
      shiftedTuple[dir] = a_tuple[dir];
    }
  shiftedTuple[a_dir] += a_shift;
  return (*this)[shiftedTuple];
}
template <class T> T& MDArray<T>::operator[](int a_linearIndex) const
{
  return m_data[a_linearIndex];
}
template <class T> T& MDArray<T>::operator[](int a_tupleIndex[DIM]) const
{
  int index = m_box.linearIndex(a_tupleIndex);
  return m_data[index];
}
template <class T> const Box& MDArray<T>::getBox() const
{
  return m_box;
}
