#include "FEPoissonOperator.H"

/// Default Constructor
FEPoissonOperator::FEPoissonOperator()
{
}

/// On construction, FEPoissonOperator fills in m_matrix using the alogirthm for L from the lecture    
FEPoissonOperator::FEPoissonOperator(const FEGrid& a_grid)
{
  int numElems;
  int matrixPos[2];
  int interSize;
  float localGradient1[DIM];
  float localGradient2[DIM];
  float gradientDot;
  float elemArea;
  Node elemNode1;
  Node elemNode2;

  m_grid = a_grid;
  numElems = m_grid.getNumElts();
  /// Get the number of interior nodes there are
  interSize = m_grid.getNumInteriorNodes();
  /// Resize the matrix
  m_matrix = SparseMatrix(interSize, interSize);

  for(int elemNum=0; elemNum<numElems; elemNum++)
    {
      for(int localNode1=0; localNode1<3; localNode1++)
	{
	  for(int localNode2=0; localNode2<3; localNode2++)
	    {
	      /// Grab the local node data corresponding to elemNum
	      elemNode1 = m_grid.getNode(elemNum, localNode1);
	      elemNode2 = m_grid.getNode(elemNum, localNode2);
	      if(elemNode1.isInterior() && elemNode2.isInterior())
		{
		  /// Compute the gradient
		  m_grid.gradient(localGradient1, elemNum, localNode1);
		  m_grid.gradient(localGradient2, elemNum, localNode2);
		  /// Compute the Area	
		  elemArea = m_grid.elementArea(elemNum);
		  /// Dot product, maybe I'll write it somewhere else later
		  gradientDot = 0;
		  for(int i=0; i<DIM; i++)
		    {
		      gradientDot+=localGradient1[i]*localGradient2[i];
		    }
		  matrixPos[0] = elemNode1.getInteriorNodeID();
		  matrixPos[1] = elemNode2.getInteriorNodeID();
		  /// Append onto the matrix
		  m_matrix[matrixPos] += elemArea*gradientDot;
		}
	    }
	}
    }
}

/// the Right Hand Side is the pseudo-code for 'b' from the lectures.                                                                                       
void FEPoissonOperator::makeRHS(vector<float> & a_rhsAtNodes, const vector<float> & a_FCentroids) const
{
  /// Get the number of elements	
  int numElems = m_grid.getNumElts();
  int globalNodeID[3];
  float elemArea;
  Element elemHold;
  Node elemNode;

  /// Resize a_rhsAtNodes to be length of number of interior nodes
  a_rhsAtNodes.resize(m_grid.getNumInteriorNodes());
  /// Resetting a_rhsAtNodes
  for(unsigned int i=0; i<a_rhsAtNodes.size(); i++)
    {
      a_rhsAtNodes[i] = 0;
    }

  /// Loop over the elements
  for(int elemNum=0; elemNum<numElems; elemNum++)
    {
      /// Calculate the area of that element
      elemArea = m_grid.elementArea(elemNum);


      elemHold = m_grid.element(elemNum);
      elemHold.vertices(globalNodeID);
      /// Loop over the nodes in the element
      for(int localNode=0; localNode<3; localNode++)
	{
	  /// Get the global node corresponding to the local node of the element
	  elemNode = m_grid.getNode(elemNum, localNode);
	  /// If it is an interior node...
	  if(elemNode.isInterior())
	    {
	      /// Then calculate the following to put at the rhs
	      //a_rhsAtNodes[elemNode.getInteriorNodeID()] += elemArea * a_FCentroids[globalNodeID[localNode]] * 1/3;
	      a_rhsAtNodes[elemNode.getInteriorNodeID()] += elemArea * a_FCentroids[elemNum] * 1/3;
	    }
	}
    }
}
      

/// return the value of m_grid.                                                                                                                             
const FEGrid& FEPoissonOperator::getFEGrid() const
{
  return m_grid;
}                                                                                                                                                             
/// return m_matrix                                                                                                                                         
const SparseMatrix& FEPoissonOperator::matrix() const
{
  return m_matrix;
}
