#ifndef _FIELDDATA_H_
#define _FIELDDATA_H_

#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <list>
#include <tr1/memory>
#include "DIM.H"
#include "Real.H"
#include "Box.H"
#include "MDArray.H"
#include "FFT1D.H"
using namespace std;
class DeltaVelocity;
class FieldData
{
public:
  FieldData();
  FieldData(std::tr1::shared_ptr<FFT1D> a_fft1dptr, int a_nghosts);
  //~FieldData();
  const MDArray<double>& operator[](int a_component) const;
  MDArray<double>& operator[](int a_component);
  void fillGhosts(MDArray<double>& array);
  void fillGhosts();
  void increment(const double& a_scalar,
                 const DeltaVelocity& a_fieldIncrement);
  void copyTo(FieldData& a_FieldData) const;
  Box m_grid;
  int m_M,m_N,m_ghosts;
  MDArray<double> m_data[DIM]; 
  std::tr1::shared_ptr<FFT1D> m_fft1dptr;
  
private:
  FieldData(const FieldData& a_FieldData);
  void operator=(const FieldData& a_FieldData);

};
#endif
