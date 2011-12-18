
#include "FEGrid.H"
#include "FEPoissonOperator.H"
#include "ReInsert.H"
#include "JacobiSolver.H"

#include <iostream>


using namespace std;

float sourceFunction(float x[DIM])
{
  float val=-.2;
  // region 1
  float Rsquared=(x[1]-9)*(x[1]-9)+x[0]*x[0];
  if(Rsquared > 25 && Rsquared < 36)
    {
      val = 1.5;
    }
  return val;
}

int main(int argc, char** argv)
{
  if(argc != 2)
    {
      cout << "this program takes one argument that is the .node and .ele ";
      cout << "file prefix"<<endl;
      return 1;
    }
  string prefix(argv[1]);
  string nodeFile = prefix+".node";
  string eleFile  = prefix+".ele";

  FEGrid grid(nodeFile, eleFile);

  FEPoissonOperator op(grid);

  const SparseMatrix& A = op.matrix(); 
  A.diagDominance();
  if(!A.symmetric())
    {
      cout<<" error in matrix assembly.  This should be a symmetric matrix\n";
      return 2;
    }
  int index[2] = {0,0};
  
  for(unsigned int i=0; i<A.M(); i++, index[0]++, index[1]++)
    {
      if(A[index] <=0)
	{
	  cout <<"negative or zero diagonal detected "<<index[0]<<endl;
	}
    }

  int nElements = grid.getNumElts();
  vector<float> sourceTerms(nElements);
  float centroid[DIM];
  for(int i=0; i<nElements; i++)
    {
      grid.centroid(centroid, i);
      sourceTerms[i] =sourceFunction(centroid);
    }

  vector<float> rhs, internalNodes, phi;
  op.makeRHS(rhs, sourceTerms);

  JacobiSolver solver;
  int iterations = 350;
  float residual = solver.solve(op.matrix(), rhs, 1E-3, iterations, internalNodes);

  reinsert(grid, internalNodes, phi);

  FEWrite(&grid, &phi, "solution.vtk");

  cout<<" Final Solver residual was "<<residual<<endl;
  
  return 0;
  
}
