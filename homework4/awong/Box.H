#ifndef _BOX_H_
#define _BOX_H_
#define DIM 2
class Box
{
public:
  Box();
  /// constructor for nontrivial Box. 
  Box(const int a_lowCorner[DIM],const int a_highCorner[DIM]);
  /// copy constructor.
  Box(const Box& a_Box);
  ///  Computes the intersection of the Box with rightBox.
  Box operator*(const Box& a_rightBox) const;
  /// Computes shifted Box by a_offset in direction a_direction.
  Box shift(int a_direction, int a_offset) const;
  /// Computes Box shifted by a_tuple.
  Box shift(int a_tuple[DIM]) const;
  /// Grow in all of the coordinate directions by a_numpoints.
  Box grow(int a_numpoints) const;
  /// Grow in coordinate direction dir by a_numpoints.
  Box grow(int dir, int a_numpoints) const;
  /// Access functions for low, high bounds of Box.
  void getLowCorner(int a_lowercorner[DIM]) const;
  void getHighCorner(int a_lowercorner[DIM]) const;
  /// returns a Box corresponding to the low (dir <0) or high (dir > 0) boundary
  /// in the |dir| coordinate direction.
  Box getBoundary(int dir) const;
  /// Computes a linear index corresponding to the input array. 
  int linearIndex(const int a_tupleIndex[DIM]) const;
  /// Computes a tuple corresponding to a linear location.
  void tupleIndex(int a_linearIndex, int a_tupleIndex[DIM]) const;
  /// Computes the size of the Box.
  int sizeOf() const;
  /// Computes length of Box on a side.
  int length(int dir) const;
  bool operator==(const Box& a_rhsBox) const;
  /// checks for empty box.
  const bool& isEmpty() const;
  /// checks for square / cubic / ...box
  bool isSquare() const;
  /// prints Box.
  void print() const;
private:  
  int m_lowCorner[DIM],m_highCorner[DIM];
  bool m_isEmpty;
};
#endif
