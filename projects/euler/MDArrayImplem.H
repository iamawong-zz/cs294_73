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
  m_data = NULL;
  m_isDefined=false;
}
template <class T> void MDArray<T>::define(Box a_box)
{
  if (m_isDefined)
    {
      delete [] m_data;
    }
  m_box = a_box;
  m_data = new T[m_box.sizeOf()];
  m_isDefined=true;
}
template <class T> void MDArray<T>::setVal(T a_val)
{
  for (int k=0;k < m_box.sizeOf();k++)
    {
      m_data[k] = a_val;
    }
}
template <class T> void MDArray<T>::copyTo(MDArray<T>& a_destArray) const
{
  if (m_box == a_destArray.m_box)
    {
      for (int k = 0;k < m_box.sizeOf();k++)
        {
          a_destArray.m_data[k] = m_data[k];
        }
    } else
        {
        Box bi = m_box*a_destArray.m_box;
        if (!(bi.isEmpty())) 
          {
            int tuple[DIM];
            for (int k = 0; k < bi.sizeOf(); k++)
              {
                bi.tupleIndex(k,tuple);
                a_destArray[tuple] = (*this)[tuple];
              }
          }
        }
}
template <class T> void MDArray<T>::shift(int a_dir,int a_len)
{
  m_box = m_box.shift(a_dir,a_len);
}
template <class T> void MDArray<T>::shift(int a_tuple[DIM])
{
  for (int dir = 0; dir < DIM; dir++)
    { 
      m_box = m_box.shift(dir,a_tuple[dir]);
    }
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
  assert(a_linearIndex > -1);
  assert(a_linearIndex < m_box.sizeOf());
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
template <class T> void MDArray<T>::operator+=(const MDArray& a_rhs)
{
  Box bi = m_box*a_rhs.m_box;
  for (int k=0;k < bi.sizeOf();k++)
    {
      int tuple[DIM];
      bi.tupleIndex(k,tuple);
        (*this)[tuple] += a_rhs[tuple];
    }
}
template <class T> void MDArray<T>::operator-=(const MDArray& a_rhs)
{
   Box bi = m_box*a_rhs.m_box;
  for (int k=0;k < bi.sizeOf();k++)
    {
      int tuple[DIM];
      bi.tupleIndex(k,tuple);
        (*this)[tuple] -= a_rhs[tuple];
    }
}
template <class T> void MDArray<T>::operator*=(const MDArray& a_rhs)
{
  Box bi = m_box*a_rhs.m_box;
  for (int k=0;k < bi.sizeOf();k++)
    {
      int tuple[DIM];
      bi.tupleIndex(k,tuple);
        (*this)[tuple] *= a_rhs[tuple];
    }
}  

template <class T> void MDArray<T>::operator/=(const MDArray& a_rhs)
{
  Box bi = m_box*a_rhs.m_box;
  for (int k=0;k < bi.sizeOf();k++)
    {
      int tuple[DIM];
      bi.tupleIndex(k,tuple);
        (*this)[tuple] /= a_rhs[tuple];
    }
} 
template <class T> void MDArray<T>::operator+=(const T& a_rhs)
{
  for (int k=0;k < m_box.sizeOf();k++)
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

