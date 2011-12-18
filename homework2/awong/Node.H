#ifndef _NODE_H_
#define _NODE_H_
#define DIM 2
#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>
using namespace std;
class Node
{
 public:
  /// default constructor. 
  Node();
  /// Real constructor: location in space and flag indicating whether node is interior
  /// or on boundary.
  Node(float a_position[DIM],
       const int&  a_interiorNodeID, 
       const bool& a_isInterior);
  /// Constant access to node Location in space.
  void  getPosition(float a_position[DIM]) const;
  const int& getInteriorNodeID() const;
  const bool& isInterior() const;

 private:
  float m_position[DIM];
  bool  m_isInterior;
  int m_interiorNodeID;
};
#endif
