#include <iostream>
#include "JacobiSolver.H"

using namespace std;

float JacobiSolver::solve(const SparseMatrix& a_A, const vector<float>& a_rhs, float a_tolerance, int a_iter, vector<float>& a_phi)
{
  /// Used to keep track of the residual is greater than tolerance
  bool greaterThanTol = true;
  /// Used to keep track of the # of iterations
  int iterNum=0;
  /// Used to find maxL[k,k]
  int indexing[2];
  float alpha = 0.85;
  float rhsNorm;
  float lhsNorm=0;
  float lambda;
  float maxDiag=0;
  float residual;
  vector<float> residualCheck(a_rhs.size());
  vector<float> product(a_rhs.size());
  vector<float> postProduct(a_rhs.size());
  
  /// Calculating and outputting the intial norm
  for(unsigned int i=0; i<a_rhs.size(); i++)
    {
      rhsNorm += (a_rhs[i]*a_rhs[i]);
    }
  rhsNorm = sqrt(rhsNorm);
  cout<<"The initial norm is "<< rhsNorm<<".\n"<<endl;

  /// Calculating lambda
  for(unsigned int i=0; i<a_A.M(); i++)
    {
      indexing[0]=i; 
      indexing[1]=i;
      /// checks to see if the diagonal is greater than the current diagonal
      if(a_A[indexing]>maxDiag)
	{
	  maxDiag = a_A[indexing];
	}
    }
  lambda = alpha/maxDiag;
  
  /// Initialize a_phi to be all zeros and set it's length to be number of interior nodes
  a_phi.resize(a_rhs.size());
  for(unsigned int i=0; i<a_phi.size(); i++)
    {
      a_phi[i]=0;
    }
  /// Loop if less then the number of iterations and greaterTol is true
  while(iterNum<a_iter && greaterThanTol)
    {
      rhsNorm = 0;
      lhsNorm = 0;
      /// Multiply L*a^l	
      product = a_A*a_phi;
      for(unsigned int i=0; i<a_phi.size(); i++)
	{
	  /// Find a^(l+1) = a^l + lambda * (b-L*a^l)
	  a_phi[i] = a_phi[i] + lambda*(a_rhs[i]-product[i]);
	}
      /// calculate L*a^(l+1)
      postProduct = a_A*a_phi;
      for(unsigned int i=0; i<a_phi.size(); i++)
	{
	  /// finding the residuals
	  residualCheck[i] = a_rhs[i] - postProduct[i];
	  lhsNorm+= residualCheck[i] * residualCheck[i];
	  rhsNorm += a_rhs[i]*a_rhs[i];
	}
      lhsNorm = sqrt(lhsNorm);
      rhsNorm = sqrt(rhsNorm);
      residual = lhsNorm/rhsNorm;
      /// Check if the residual is less than tolerance, if it is then return the residual otherwise do another iteration
      if(residual <= a_tolerance)
	{
	  greaterThanTol = 0;
	  cout<<"The final norm is "<< residual<<".\n"<<endl;
	  cout<<"It took "<<iterNum<<" iterations.\n"<<endl;
	  return residual;
	}
      else
	{
	  iterNum++;
	}
    }
}

