
#ifndef _REINSERT_H_
#define _REINSERT_H_
#include "FEGrid.H"
#include <vector>

/** given a solution defined on the internal nodes of the grid, 
    builds a full nodal solution with a value of 0 for the boundary nodes.
    a_grid and a_internalSolution are the inputs, a_fullVector is the output
*/
void reinsert(const FEGrid& a_grid, const vector<float>& a_internalSolution, vector<float>& a_fullVector);

#endif
