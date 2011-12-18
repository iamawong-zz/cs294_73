#ifndef _ELEMENT_H_
#define _ELEMENT_H_
#define VERTICES 3

class Element
{
public: 
  Element();
  /// Constructor. 
  Element(int a_vertices[VERTICES]);

  /// local indexing to get nodeNumber.
  const int& operator[](const int& a_localNodeNumber) const;

  void vertices(int a_vertices[VERTICES]) const;

private:
  int m_vertices[VERTICES];
};
#endif
