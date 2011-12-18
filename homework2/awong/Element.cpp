#include <cassert>
#include "Element.H"

Element::Element(int a_vertices[VERTICES])
{
  for (int ivert = 0;ivert < VERTICES;ivert++)
    {
      m_vertices[ivert] = a_vertices[ivert];
    }
};


Element::Element()
{
  for (int ivert = 0;ivert < VERTICES;ivert++)
    {
      m_vertices[ivert] = -1;
    }
};
const int& Element::operator[](const int& a_localNodeNumber) const
{
  assert(a_localNodeNumber < VERTICES);
  return m_vertices[a_localNodeNumber];
};


void Element::vertices(int a_vertices[VERTICES]) const
{
  for(int i=0; i<VERTICES; i++)
    {
      a_vertices[i] = m_vertices[i];
    }
}
