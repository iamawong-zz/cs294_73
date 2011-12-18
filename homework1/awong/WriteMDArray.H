
#ifndef _WRITEMDARRAY_H_
#define _WRITEMDARRAY_H_
#include "MDArray.H"
#include <cstring>

const char* MDWrite(MDArray* a_array);

void MDWrite(const char* filename, MDArray* a_array);

#endif 
