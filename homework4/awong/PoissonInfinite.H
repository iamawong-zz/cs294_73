#ifndef _POISSONINFINITE_H_
#define _POISSONINFINITE_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "MDArray.H"
#include "FFTW1DSine.H"
#include "FFTMDSine.H"
#define STENCILSIZE 9
using namespace std;
class PoissonInfinite
{
public:
  PoissonInfinite();
  /// solver for [0 ... a_N]^DIM grid, with infinite-domain BCs. 
  /// a_N = P*2^M M>=3, P highly composite, a_growthfactor = a_N/8.
  PoissonInfinite(unsigned int a_N, double a_h, unsigned int a_growthFactor);
  void solve(MDArray<double>& a_solution, 
        const MDArray<double>& a_rhs) const;
  /// You will need applyOperator in MLC calculation. 
  void applyOperator(MDArray<double>& a_LOfPhi, 
                     const MDArray<double>& a_phi,
                     const Box& a_bxApply) const;

protected:
  void dirichletSolve(MDArray<double>& a_solution,bool a_isInner) const;
  void computeInnerBoundaryCharges(MDArray<double> a_charges[DIM][2],
                                   const MDArray<double>& a_field) const;
  void 
  computeOuterBoundaryCharges(MDArray<double> a_outerCharges[DIM][2],
                       const MDArray<double> a_innerCharges[DIM][2]) const;
  void
  faceToFace(MDArray<double>& a_bigRHS,
             const MDArray<double>& a_array,
             const Box& a_outerBoundary) const;
  void getStencil(int a_stencil[STENCILSIZE][DIM],int tuple[DIM]) const;
  int m_N;
  double m_h,m_L;
  FFTMDSine* m_fftmd;
  FFTMDSine* m_fftmdBig;
  Box m_box,m_boxBig;
  int m_weights[STENCILSIZE];
  int m_offsets[STENCILSIZE][DIM];
};
#endif
