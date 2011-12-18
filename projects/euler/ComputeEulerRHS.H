#ifndef _COMPUTE_EULER_RHS_H_
#define _COMPUTE_EULER_RHS_H_
#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "MDArray.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include "FieldData.H"
#include "DeltaVelocity.H"
#include "Projection.H"
#include "AdvectionOperator.H"

class ComputeEulerRHS

{
public:
	void operator() (
		DeltaVelocity& a_newDv,
			const double& a_time, const double& a_dt,
			const FieldData& a_velocity,
			DeltaVelocity& a_oldDv);
};
#endif

