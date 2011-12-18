/*
 * Author: Johnny Lee
 * SVN Login: itzfx
 */

#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "Real.H"
#include "MDArray.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include "FieldData.H"
#include "DeltaVelocity.H"
#include "PoissonSolver.H"
#include "DIM.H"
#include "Real.H"
#include "Projection.H"
#include "WriteMDArray.H"

using namespace std;

Projection::Projection(std::tr1::shared_ptr<FFT1D> a_fft1dPtr)
{
  // create the solver
  m_solver = PoissonSolver(a_fft1dPtr);
  // store N
  m_N = a_fft1dPtr->getN();
}

void Projection::applyProjection(FieldData& a_velocity) 
{
 // this is const so it should only modify a_velocity
  /*
   * From the project handout, we have:
   * P = I - grad(DELTA^-1)div
   *
   * Discretized, we have:
   * P^h = I - (grad^h)((DELTA^h)^-1)div^h
   * (DELTA^h)(phi^h)_i = (1/(4h^2))(-4phi_i + phi_(i+e^0) + phi_(i-e^0) + phi_(i+e^1) + phi_(i-e^1))
   *
   * After some conversation, it seems we will use:
   * Box vBox(a_velocity.m_grid)
   * MDArray<Real> R(vBox)               <---- prepare for divergence calculation
   * divergence(R, a_velocity)           <---- stores divergence in R
   * m_solver.solve(R)                   <---- gives (DELTA^-1 div) in R
   * MDArray<Real> Rghost(vBox.grow(1))  <---- prepare for gradient calculation
   * R.copyTo(Rghost)                    <---- set the data in Rghost
   * DeltaVelocity dv;                   <---- prepare for gradient
   * dv.init(a_velocity)                 <---- prepare for gradient
   * gradient(dv, Rghost)                <---- get the gradient
   * 
   * Now we have grad(DELTA^-1)div
   * We need to subtract this from the identity matrix.
   * Another way to put this is that we need to go through each index in dv and 
   * perform 1-dv_element if the index is on the standard diagonal.
   * 
   * This is slightly tricky since the indexing is a bit strange.
   */

  // create a box that's the same as the one in a_velocity for use
  Box vBox(a_velocity.m_grid);
  a_velocity.fillGhosts();

  // Create an MDArray for div, and get the divergence
  MDArray<Real> R(vBox);
  divergence(R, a_velocity);

  // Get DELTA^-1 div using the PoissonSolver
  m_solver.solve(R);

  // now create a ghosted box for the gradient
  MDArray<Real> Rghost(vBox.grow(1));
  // and copy over DELTA^-1 div
  R.copyTo(Rghost);
  // need to fill ghosts on Rghost
  a_velocity.fillGhosts(Rghost);
  // And finally, get the gradient
  DeltaVelocity dv;
  dv.init(a_velocity);
  gradient(dv, Rghost);
  //a_velocity[0] -= dv[0];
  //a_velocity[1] -= dv[1];
  a_velocity.increment(-1.0, dv);
  //  MDWrite(dv[0]);
  
}

void Projection::gradient(DeltaVelocity& a_vector, const MDArray<Real>& a_scalar) const
{
  // http://en.wikipedia.org/wiki/Gradient
  /*
   * From the project handout, we have:
   * grad(phi) = (delta(p)/delta(x_0), ..., delta(p)/delta(x_(D-1)))
   * Discretized, we have:
   * grad^h(phi^h)_i = ((1/2h)(phi_(i+e^0) - phi_(i-e^0)), (phi_(i+e^1) - phi_(i-e^0)))
   */

  // attempt:
  // the low corner of the box
  int lowCorner[DIM];
  a_vector.m_grid.getLowCorner(lowCorner);
  // the high corner of the box
  int highCorner[DIM];
  a_vector.m_grid.getHighCorner(highCorner);
  // this is mesh spacing
  double h = 1.0 / m_N;
  Box bx = a_vector.m_grid;

  // TODO is this right???
  // I'm going to assume DIM = 2 since we're working in 2-d space
  for (int k = 0; k < bx.sizeOf();k++)
  {
    int index[2];
    bx.tupleIndex(k,index);
    int i = index[0];
    int j = index[1];
    int plusE0[2] = {i+1, j};  // i+e^0 (right side)
    int minusE0[2] = {i-1, j}; // i-e^0 (left side)
    int plusE1[2] = {i, j+1};  // i+e^1 (top side)
    int minusE1[2] = {i, j-1}; // i-e^1 (bottom side)
    
    // figure out the x portion and save it
    (a_vector[0])[index] = (a_scalar[plusE0] - a_scalar[minusE0])/(2 * h);
    // and the y portion and save it
    (a_vector[1])[index] = (a_scalar[plusE1] - a_scalar[minusE1])/(2 * h);
  }
}

void Projection::divergence(MDArray<Real>& a_scalar, const FieldData& a_vector) const
{
  // http://en.wikipedia.org/wiki/Divergence
  /*
   * From the project handout, we have:
   * div(w) = sum(d=0 to D-1, delta(w_d)/delta(x_d))
   * Discretized, we have:
   * div^h(w^h)_i = (1/2h)((u_0)_(i+e^0) - (u_0)_(i-e^0)) + (1/2h)((u_1)_(i+e^1) - (u_1)_(i-e^1))
   * Note:
   *   w^h = {u_0,u_1}, e^0 = {1,0} , e^1 = {0,1}
   */

  // ghosts have already been filled
  // the low corner of the box
  Box bx = a_vector.m_grid;
 
  // fix the corners
  // this is mesh spacing
  double h = 1.0 / a_vector.m_N;

  // TODO is this right???
  // go through every element in a_vector which should be w^h
  // I'm going to assume DIM = 2 since we're working in 2-d space
  for (int k = 0; k < bx.sizeOf();k++)
  {
    int index[2];
    bx.tupleIndex(k,index);
    int i = index[0];
    int j = index[1];
      int plusE0[2] = {i+1, j};  // i+e^0 (right side)
      int minusE0[2] = {i-1, j}; // i-e^0 (left side)
      int plusE1[2] = {i, j+1};  // i+e^1 (top side)
      int minusE1[2] = {i, j-1}; // i-e^1 (bottom side)

      // figure out the x portion
      Real xPortion = a_vector.m_data[0][plusE0] - a_vector.m_data[0][minusE0];
      // and the y portion
      Real yPortion = a_vector.m_data[1][plusE1] - a_vector.m_data[1][minusE1];
      // save the result which is 1/2h(xPortion + yPortion)
      a_scalar[index] =  (xPortion + yPortion)/(2 * h);
    
  }
}
