#ifndef _SPARSEMATRIX_H_
#define _SPARSEMATRIX_H_

#include <vector>
#include <cassert>
#include <cmath>

using namespace std;
class SparseMatrix
{
public:
 
  /// a matrix with zero rows and zero columns
  SparseMatrix();

  /// set up an M rows and N columns sparse matrix with all values of zero (no non-zero elements)
  // this should make m_data and m_colIndex the correct size, but not create any matrix fill yet.
  // don't forget to initialize m_zero 
  SparseMatrix(unsigned int a_M, unsigned int a_N);

  /// Matrix Vector multiply.  a_v.size()==N, returns vector of size M
  vector<float> operator*(const vector<float>& a_v) const;

  ///accessor functions for get and set operations of matrix elements
  // If the element does not already exist in the matrix create storage for it and return a reference
  float& operator[](int a_index[2]);

  ///accessor function just to get a value
  // if the element does not already exist in the matrix, then return a reference to m_zero
  const float& operator[](int a_index[2]) const;

  /// zero out all the elements, but leave the sparse structure in place.
  // this modified m_data to hold the numeric value of zero.  m_colIndex is unchanged
  void zero();

  /// checks if the matrix is diagonally dominated
  const  void diagDominance() const;

  /// build and return a new SparseMatrix that is the transpose of the input matrix.
  SparseMatrix transpose() const;

  /// return number of rows
  unsigned int M() const;

  /// return number of columns
  unsigned int N() const;

  /// Check that this matrix is symmetrix A(i,j)==A(j,i) for all i and j in matrix.
  // careful not to accidentally create a dense matrix with your accessors.
  bool symmetric() const;

  /// dump something human readable to the screen.  Your choice how this looks
  void print() const;

private:
  unsigned int m_m, m_n;
  float m_zero;
  vector<vector<float> > m_data;
  vector<vector<int> >   m_colIndex;
};

#endif
