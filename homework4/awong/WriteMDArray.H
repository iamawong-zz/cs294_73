
#ifndef _WRITEMDARRAY_H_
#define _WRITEMDARRAY_H_
#include "MDArray.H"
class ParticleSet;

#include <cstring>

const char* MDWrite(MDArray<float>* a_array);

const char* MDWrite(MDArray<double>& a_array);

void MDWrite(const char* filename, MDArray<float>* a_array);

void MDWrite(const char* filename, MDArray<double>& a_array);

void PWrite(const char* filename, const ParticleSet* a_particles);

const char* PWrite(const ParticleSet* a_particles);

#endif 
