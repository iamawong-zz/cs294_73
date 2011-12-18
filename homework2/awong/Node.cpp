#include <cstdio>
#include <cmath>
#include <cassert>
#include <cfloat>
#include <vector>
#include "Node.H"
using namespace std;
Node::Node()
{
  for (int idir = 0;idir < DIM;idir++)
    {
      m_position[idir] = FLT_MAX;
    }
  m_isInterior = true;
  m_interiorNodeID = -1;
};
  
Node::Node(float a_position[DIM],
           const int& a_interiorNodeID,
           const bool& a_isInterior)
{
  for (int idir = 0;idir < DIM;idir++)
    {
      m_position[idir] = a_position[idir];
    }
  m_isInterior = a_isInterior;
  m_interiorNodeID = a_interiorNodeID;
};

void Node::getPosition(float a_x[DIM]) const
{
  for (int idir=0;idir < DIM; idir++)
    {
      a_x[idir] = m_position[idir];
    }
};

const int& Node::getInteriorNodeID() const
{
  return m_interiorNodeID;
};

const bool& Node::isInterior() const
{
  return m_isInterior;
};


  
