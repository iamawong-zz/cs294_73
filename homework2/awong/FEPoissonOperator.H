#ifndef _FEPOISSONOPERATOR_H_
#define _FEPOISSONOPERATOR_H_
#include <cstdio> 
#include <cmath> 
#include <cassert>
#include <cmath>
#include <vector>
using namespace std;

#include "Node.H"
#include "Element.H"
#include "FEGrid.H"
#include "SparseMatrix.H"

/// A overal Class comment can go here.  FEPoissonOperator creates the L matrix and b vector from class
class FEPoissonOperator
{
 public:
  /// null constructor, empty SparseMatrix, empty FEGrid
  FEPoissonOperator();

  /// On construction, FEPoissonOperator fills in m_matrix using the alogirthm for L from the lecture
  FEPoissonOperator(const FEGrid& a_grid);

  /// the Right Hand Side is the pseudo-code for 'b' from the lectures.
  void makeRHS(vector<float> & a_rhsAtNodes, const vector<float> & a_FCentroids) const;

  /// return the value of m_grid.  
  const FEGrid& getFEGrid() const;

  /// return m_matrix
  const SparseMatrix& matrix() const;

 private:
  SparseMatrix m_matrix;
  FEGrid m_grid;
};
#endif 
