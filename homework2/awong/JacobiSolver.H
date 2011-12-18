#ifndef _JACOBISOLVER_H_
#define _JACOBISOLVER_H_

#include "SparseMatrix.H"
#include <vector>

class JacobiSolver
{
public:
  /** solves until max norm of residual is less than tolerance, 
      or until a_iter iterations have been completed, whichever happens first.  
      returns final residual ratio in max norm ||b-Ax||/||b||
      the solution is stored back into a_phi.
      a_A is dimensioned as numInterNodes, as is a_rhs, and a_phi
  */
  float solve(const SparseMatrix& a_A, const vector<float>& a_rhs, float a_tolerance, int a_iter, vector<float>& a_phi);
};

#endif


