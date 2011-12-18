#ifndef _FEGRID_H_
#define _FEGRID_H_

#include <cstdio>   
#include <cmath>    
#include <cassert>
#include <cstring> 
#include <vector> 
#include <string>
#include "Node.H"   
#include "Element.H"

using namespace std;

class FEGrid
{
public:
  FEGrid();
  /// Constructor by reading from file.
  FEGrid(const std::string& nodeFile, const std::string& a_elementFileName);

  /// Compute gradient of basis function at node a_nodeNumber = 0,..,DIM, on element a_eltNumber. 
  void gradient(float a_gradient[DIM], 
		const int& a_eltNumber,
		const int& a_localNodeNumber) const;
  /// Compute centroid of element.                                
  void centroid(float a_centroid[DIM],
		const int& a_eltNumber) const;
  /// Compute area of element.
  float elementArea(const int& a_eltNumber) const;
 
  /// get reference to node                      
  const Node& getNode(const int& a_eltNumber,
                      const int& a_localNodeNumber) const; 
  /// Get number of elements.
  int getNumElts() const;
  int getNumNodes() const;
  int getNumInteriorNodes() const;

  const Element& element(int i) const;
  const Node& node(int i) const;
  private:                            
  vector<Node > m_nodes;
  vector<Element > m_elements;
  int m_numInteriorNodes;
};

const char* FEWrite(FEGrid* a_grid, const char* a_filename);
const char* FEWrite(FEGrid* a_grid);
const char* FEWrite(FEGrid* a_grid, vector<float>* a_scalarField, const char* a_filename);
const char* FEWrite(FEGrid* a_grid, vector<float>* a_scalarField);


#endif
