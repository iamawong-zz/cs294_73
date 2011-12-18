#ifndef _ADVECTIONOPERATOR_H_
#define _ADVECTIONOPERATOR_H_
#include <iostream>
#include <cmath>
#include <vector>
#include "Box.H"
#include "MDArray.H"
#include "FieldData.H"
#include "DeltaVelocity.H"
void advectionOperator(DeltaVelocity& a_divuu,
                       const FieldData& a_velocity,
                       const Box m_grid,
                       const double& a_h);
#endif
