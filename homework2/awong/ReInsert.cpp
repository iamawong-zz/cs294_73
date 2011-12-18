#include "ReInsert.H"

void reinsert(const FEGrid& a_grid, const vector<float>& a_internalSolution, vector<float>& a_fullVector)
{
  /// Make sure that the size of the full solution is equal to the number of total nodes
  assert(a_internalSolution.size() == a_grid.getNumInteriorNodes());
  /// Resize the full vector solution
  a_fullVector.resize(a_grid.getNumNodes());

  Node nodeLooper;
  /// Loop over all the nodes
  for(int i=0; i<a_grid.getNumNodes(); i++)
    {
      /// nodeLooper is to grab each node
      nodeLooper = a_grid.node(i);
      if(nodeLooper.isInterior())
	{
	  /// If the node is an interior node, then we copy the internal solution over
	  a_fullVector[i] = a_internalSolution[nodeLooper.getInteriorNodeID()];
	}
      else
	{
	  /// If not, we set it to zero
	  a_fullVector[i] = 0;
	}
    }
}
      
      
