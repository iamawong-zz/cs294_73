#include <vector>
#include "Box.cpp"

using namespace std;

void square_dgemm(int n, double *A, double *B, double *C)
{
  double A_i[n];
  double B_j[n];
  double C_ij; 
  
  for(int i=0; i<n; i++)
    {
      /// read row i into A_i for fast memory
      for(int m=0; m<n; m++)
	{
	  A_i[m] = A[i+(m*n)];
	}
      for(int j=0; j<n; j++)
	{
	  /// Read C_ij into memory slot
	  C_ij = C[(n*j)+i];
	  /// read column j into B_j for fast memory
	  for(int m=0; m<n; m++)
	    {
	      B_j[m] = B[(j*n)+m];
	    }
	  for(int k=0; k<n; k++)
	    {
	      C_ij = C_ij + (A_i[k] * B_j[k]);
	    }
	  /// Rewrite back into 'slow memory'
	  C[(n*j)+i] = C_ij;
	}
    }
}
