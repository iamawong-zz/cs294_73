#include "AdvectionOperator.H"

void advectionOperator(DeltaVelocity& a_divuu, 
		       const FieldData& a_velocity,
		       const Box a_grid,
		       const double& a_h)
{
  int accessTuple[DIM];
  int shiftDownTuple[DIM];
  int shiftUpTuple[DIM];
  FieldData currVelocity;
  a_velocity.copyTo(currVelocity);
  currVelocity.fillGhosts();

  /// This is to traverse the grid
  //// INITIALIZING parameters a_divuu
  a_divuu.init(a_velocity);
  for(int i=0; i<a_grid.length(0); i++)
    {
      accessTuple[0]=i;
      for(int j=0; j<a_grid.length(1); j++)
	{
	  accessTuple[1]=j;
	  /// k is to go through the components a_divuu_k
	  for(int k=0; k<DIM; k++) {
	    for(int l=0; l<DIM; l++) {
	      double left;
	      double right;
	      double deriv;
	      
	      shiftUpTuple[0] = accessTuple[0];
	      shiftUpTuple[1] = accessTuple[1];
	      shiftDownTuple[0] = accessTuple[0];
	      shiftDownTuple[1] = accessTuple[1];
	      /// Used for shifting to the index on the left or right
	      shiftUpTuple[l]+=1;
	      shiftDownTuple[l]-=1;
	      
	      left = (currVelocity.m_data[k][accessTuple]+currVelocity.m_data[k][shiftUpTuple])*(a_velocity.m_data[l][accessTuple]+a_velocity.m_data[l][shiftUpTuple])/4;
	      right = (currVelocity.m_data[k][accessTuple]+currVelocity.m_data[k][shiftDownTuple])*(currVelocity.m_data[l][accessTuple]+currVelocity.m_data[l][shiftDownTuple])/4;
	      
	      deriv = (left-right)/a_h;
	      // adds onto divuu_k
	      a_divuu.m_data[k][accessTuple] += deriv;
	    }
	  }
	}
    }
}
