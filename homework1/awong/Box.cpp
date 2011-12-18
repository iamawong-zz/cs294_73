#include "Box.H"
#include <cmath>
#include <cassert>
#include <cstdio>

/// default constructor
Box::Box()
{}

/// constructor with given dimensions
Box::Box(const int a_lowCorner[DIM], const int a_highCorner[DIM])
{
  for(int i=0; i<DIM;i++)
    {
      m_lowCorner[i] = a_lowCorner[i];
      m_highCorner[i] = a_highCorner[i];
    }
}

/// copy constructor
Box::Box(const Box& a_Box)
{
  a_Box.getLowCorner(m_lowCorner); /// grab the coordinates
  a_Box.getHighCorner(m_highCorner);
}

/// intersection of two boxes
Box Box::operator*(const Box& a_rightBox) const
{
  int rhs_lowCorner[DIM];
  int rhs_highCorner[DIM];
  int intersect_lowCorner[DIM];
  int intersect_highCorner[DIM];

  a_rightBox.getLowCorner(rhs_lowCorner);
  a_rightBox.getHighCorner(rhs_highCorner);
  for(int i=0; i<DIM; i++)
    {
      /// comparison of low corner of box B with the corners to see if its in between
      if(rhs_lowCorner[i]>m_lowCorner[i] && rhs_lowCorner[i]<m_highCorner[i])
	{
	  intersect_lowCorner[i] = rhs_lowCorner[i];
	}
      /// comparison of low corner of box A with the corners to see if its in between
      else if(m_lowCorner[i]>rhs_lowCorner[i] && m_lowCorner[i]<rhs_highCorner[i])
	{
	  intersect_lowCorner[i] = m_lowCorner[i];
	}
      /// If no intersection, then okay
      else
	{
	  /// stop the program if no intersection
	  assert(0);
	  intersect_lowCorner[i]=0;
	}
      /// comparison of the high corner of box B with the corners to see if it is in between
      if(rhs_highCorner[i]>m_lowCorner[i] && rhs_highCorner[i]<m_highCorner[i])
	{
	  intersect_highCorner[i] = rhs_highCorner[i];
	}
      /// comparison of the high corner of box A with the corners to see if it is in between
      else if(m_highCorner[i]>rhs_lowCorner[i] && m_highCorner[i]<rhs_highCorner[i])
	{
	  intersect_highCorner[i]=m_highCorner[i];
	}
      else
	{
	  /// stop the program if no intersection
	  assert(0);
	  intersect_highCorner[i]=0;
	}
    }
    
  return Box(intersect_lowCorner, intersect_highCorner);
}

/// shifting the box in a_direction by a_offset
Box Box::shift(int a_direction, int a_offset)
{
  assert(a_direction>0 && a_direction<DIM);
  int newlowCorner[DIM];
  int newhighCorner[DIM];

  /// Copy them over, and if it's the right dimension then add the offset
  for(int i=0;i<DIM;i++)
    {
      if(i==a_direction)
	{
	  newlowCorner[i] = m_lowCorner[i]+a_offset;
	  newhighCorner[i] = m_highCorner[i]+a_offset;
	}
      else
	{
	  newlowCorner[i] = m_lowCorner[i];
	  newhighCorner[i] = m_highCorner[i];
	}
    }

  return Box(newlowCorner, newhighCorner);
}

/// grow a shell around your box by a_numpoints
Box Box::grow(int a_numpoints)
{
  int newlowCorner[DIM];
  int newhighCorner[DIM];

  /// loop through the points, and do the operation on them
  for(int i=0; i<DIM; i++)
    {
      newlowCorner[i] = m_lowCorner[i]-a_numpoints;
      newhighCorner[i] = m_highCorner[i]+a_numpoints;
    }

  return Box(newlowCorner, newhighCorner);
}

/// get the lower corner and return the input tuple
void Box::getLowCorner(int a_lowercorner[DIM]) const
{

  for(int i=0; i<DIM; i++)
    {
      a_lowercorner[i] = m_lowCorner[i];
    }
  
  return;
}

/// get higher corner and return the input tuple
void Box::getHighCorner(int a_higherCorner[DIM]) const
{
  for(int i=0; i<DIM; i++)
    {
      a_higherCorner[i] = m_highCorner[i];
    }
  
  return;
}

/// gives a linear index based on the tuple that we want, implementation only for DIM=2
int Box::linearIndex(const int a_tupleIndex[DIM]) const
{
  return a_tupleIndex[1] - m_lowCorner[0] + (a_tupleIndex[0] - m_lowCorner[1])*(m_highCorner[0] - m_lowCorner[0] + 1);
}

/// given a linear index this finds the tuple location and returns it in a_tupleIndex
/// currently only works for 2D, and it is right
void Box::tupleIndex(int a_linearIndex, int a_tupleIndex[DIM]) const
{
  a_tupleIndex[1] = a_linearIndex%(m_highCorner[0]-m_lowCorner[0]+1)+m_lowCorner[0];
  a_tupleIndex[0] = a_linearIndex/(m_highCorner[0]-m_lowCorner[0]+1)+m_lowCorner[1];
}

/// finds the area/volume/hypervolume of box
int Box::sizeOf() const
{
  int dimensions[DIM];
  int size=1;

  for(int i=0; i<DIM; i++)
    {
      /// dimensions calculates the length of each dimension
      dimensions[i] = m_highCorner[i] - m_lowCorner[i] + 1;
      /// multiply the lengths together
      size*=dimensions[i];
    }

  return size;
}

/// comparison of box A and box B
bool Box::operator==(const Box& a_rhsBox) const
{
  bool equals=1;
  int rhs_lowCorner[DIM];
  int rhs_highCorner[DIM];

  a_rhsBox.getLowCorner(rhs_lowCorner);
  a_rhsBox.getHighCorner(rhs_highCorner);

  for(int i=0; i<DIM; i++)
    {
      /// if the corners are matching, keep going
      while(equals)
	{
	  /// checks if the lower corner and higher corner are matching
	  if(m_lowCorner[i]==rhs_lowCorner[i] && m_highCorner[i]==rhs_highCorner[i])
	    {}
	  /// Or else stop the loop
	  else
	    {
	      /// stops the while loop
	      equals=0;
	    }
	}
    }

  return equals;
}
