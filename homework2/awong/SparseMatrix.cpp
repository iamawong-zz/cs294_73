#include <iostream>
#include "SparseMatrix.H"

using namespace std;
/// ostream overload for a vector of template
/*
template <class T>
ostream& operator<<(ostream& a_os, const vector<T>& a_vec);
*/
/// ostream overload for a vector of floats
ostream& operator<<(ostream& a_os, const vector<float>& a_vec);
/// ostream overload for a vector of ints
ostream& operator<<(ostream& a_os, const vector<int>& a_vec);

ostream& operator<<(ostream& a_os, const vector<float>& a_vec)
{
  a_os << "[";
  for(int i=0; i<a_vec.size(); i++)
    {
      a_os << a_vec[i] << " ";
    }
  a_os << "]\n";
  return a_os;
}

ostream& operator<<(ostream& a_os, const vector<int>& a_vec)
{
  a_os << "[";
  for(int i=0; i<a_vec.size(); i++)
    {
      a_os << a_vec[i] << " ";
    }
  a_os << "]\n";
  return a_os;
}

/*ostream& operator<<(ostream& a_os, const vector<T>& a_vec)
{
  for(int i=0; i<a_vec.size(); i++)
    {
      a_os << "[" << a_vec[i] << " ";
    }
  a_os << "]\n";
  return a_os;
  }
*/

/// Default Constructor
SparseMatrix::SparseMatrix()
{
  m_zero = 0;
}

/// Constructor
SparseMatrix::SparseMatrix(unsigned int a_M, unsigned int a_N)
{
  /// Make sure it's positive numbers size
  assert(a_M > 0 && a_N > 0);
  /// Set the sizes to the input
  m_m = a_M;
  m_n = a_N;
  /// Initialize the number of rows for m_data and m_colIndex
  m_data.resize(a_M);
  m_colIndex.resize(a_M);
  /// Initialize the return m_zero 
  m_zero = 0;
}

/// Matrix-Vector multiplication
vector<float> SparseMatrix::operator*(const vector<float>& a_v) const
{
  /// Initialize a vector of floats of length m_m and set the values to zero
  vector<float> retFloats(m_m,0);
  int colIndex;

  /// Loop through the rows
  for(unsigned int i=0; i<m_colIndex.size(); i++)
    {
      /// Loop through the columns up to the length of the vector in the i-th row
      for(unsigned int j=0; j<m_colIndex[i].size(); j++)
	{
	  /// Copy the index that is needed over to the placeholder
	  colIndex = m_colIndex[i][j];
	  /// Do the multiplication and addition to the return vector
	  retFloats[i] += m_data[i][j]*a_v[colIndex];
	}
    }
  
  return retFloats;
}

/// Grab a data point with the intention to use
float& SparseMatrix::operator[](int a_index[2])
{
  /// Loop through the size of m_colIndex to see if the column is already set
  for(unsigned int j=0; j<m_colIndex[a_index[0]].size(); j++)
    {
      /// If the column is already set, return the data
      if(m_colIndex[a_index[0]][j] == a_index[1])
	{
	  return m_data[a_index[0]][j];
	}
    }
  
  assert(m_colIndex[a_index[0]].size() < this->N());
  /// If the column isn't set, push_back on the vector and add new data
  m_colIndex[a_index[0]].push_back(a_index[1]);
  m_data[a_index[0]].push_back(m_zero);
  return m_data[a_index[0]][m_colIndex[a_index[0]].size()-1];
}

const float& SparseMatrix::operator[](int a_index[2]) const
{
  /// Loop through the size of m_colIndex to see if the column is already set
  for(unsigned int j=0; j<m_colIndex[a_index[0]].size(); j++)
    {
      /// If the column is already set, return the data, otherwise return zero
      if(m_colIndex[a_index[0]][j] == a_index[1])
	{
	  return m_data[a_index[0]][j];
	}
    }
  return m_zero; 
}
  

/// Zeros out all the data, but retain the sparse structure
void SparseMatrix::zero()
{
  /// Loop through the rows
  for(unsigned int i=0; i<m_m; i++)
    {
      for(unsigned int j=0; j<m_data[i].size(); j++)
	{
	  /// Loop through the columns and set the data to zero
	  m_data[i][j]=m_zero;
	}
    }
}

/// Returns the number of rows
unsigned int SparseMatrix::M() const
{
  return m_m;
}

/// Returns the number of columns
unsigned int SparseMatrix::N() const
{
  return m_n;
}


/// Construct the transpose of the matrix
SparseMatrix SparseMatrix::transpose() const
{
  int dataPoint[2];
  /// New transpose matrix NxM
  SparseMatrix transMatrix(this->N(), this->M());
  for(int i=0; i<m_m; i++)
    {
      for(int j=0; j<m_colIndex[i].size();j++)
	{
	  /// Puts the data in (i,j) to the transpose in location (colIndex,i)
	  dataPoint[0] = m_colIndex[i][j];
	  dataPoint[1] = i;
	  transMatrix[dataPoint] = m_data[i][j];
	}
    }

  return transMatrix;
}

/// Checks if the matrix is symmetric
bool SparseMatrix::symmetric() const
{
  /// Stores the indices
  int indexEle[2];
  int indexSym[2];
  
  for(unsigned int i=0; i<m_m; i++)
    {
      /// Loop through the rows
      for(unsigned int j=0; j<m_colIndex[i].size(); j++)
	{
	  indexEle[0]=i;
	  indexEle[1]=j;
	  indexSym[0]=j;
	  indexSym[1]=i;
	  if((*this)[indexEle] != (*this)[indexSym])
	    {
	      return false;
	    }
	}
    }
  return true;
}

void SparseMatrix::print() const
{
  cout<<"This is the colIndex\n"<<endl;
  for(unsigned int i=0; i<this->M(); i++)
    {
      cout<<m_colIndex[i]<<endl;
    }
  
  cout << "This is the member data\n" <<endl;
  for(unsigned int i=0; i<this->M(); i++)
    {
      cout<<m_data[i]<<"\n"<<endl;
    }
}

/// checks if the matrix is diagonally dominated
const void SparseMatrix::diagDominance() const
{
  int diagInd[2];
  int offDiagInd[2];
  float diag;
  float offdiagSum;
  for(unsigned int i=0; i<m_m; i++)
    {
      diagInd[0]=i;
      offDiagInd[0]=i;
      diagInd[1]=i;
      diag = abs((*this)[diagInd]);
      offdiagSum = 0;
      for(unsigned int j=0; j<m_colIndex[i].size(); j++)
	{
	  /// if it is (i,i) then add set it to diag
	  if(i!=j)
	    {
	      offDiagInd[1]=j;
	      offdiagSum += abs((*this)[offDiagInd]);
	    }
	}
      if(diag<offdiagSum)
	{
	  cout<<"It is not diagonally dominant for row "<<i<<"."<<endl;
	}
    }
  cout<<"The matrix is diagonally dominant!"<<endl;
}
