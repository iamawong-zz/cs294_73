#include "MDArray.H"
#include "WriteMDArray.H"
#include <cstdio>
#include <cassert>
#include <cmath>
int main(int argc, char* argv[])
{
  /* Test program: compute the Laplacian on a grid. Grid size is input on the command line 
     as DIM integers > 2. Length of the domain in the first dimesnion is 1. */

  assert(argc == DIM+1);
  int hi[DIM];
  int lo[DIM];

  for (int i = 0; i < DIM; i++)
    {
      sscanf(argv[i+1],"%d",&hi[i]);
      lo[i] = 0;
    }
  float h = 1./(hi[0]);
  Box domainWithGhost(lo,hi);
  Box domain = domainWithGhost.grow(-1);
  MDArray phi(domainWithGhost);
  MDArray LOfPhi(domain);
  for (int i=0;i<domainWithGhost.sizeOf();i++)
    {
      int tuple[DIM];
      domainWithGhost.tupleIndex(i,tuple);
      float val = 1.;
      for (int dir = 0;dir < DIM; dir++)
        {
          if (dir%2 == 0)
            {
              val = val*sin(2*M_PI*tuple[dir]*h);
            }
          else
            {
              val = val*sin(2*M_PI*tuple[dir]*h);
            }
        } 
      phi[i] = val;
    }
  
  for (int i=0;i < domain.sizeOf();i++)
    {
      int tuple[DIM];
      LOfPhi[i] = 0.;
      domain.tupleIndex(i,tuple);
      assert (i==domain.linearIndex(tuple));
      int iphi=domainWithGhost.linearIndex(tuple);
      for (int dir = 0;dir < DIM; dir++)
        {
          LOfPhi[i] += (phi.indexShift(tuple,dir,1) + phi.indexShift(tuple,dir,-1));          
        }
       LOfPhi[tuple] -= 2*DIM*phi[tuple];
    }
  LOfPhi /= h*h;
  /* debug test: if the grid is square / cubic, LOFPhi is proportional to
     phi. We use a trick to find the constant of proportionality, then
     check to see if phi is proportional to LOFPhi to within roundoff. */
  MDWrite("phi", &phi);
  MDWrite("LOfPhi", &LOfPhi);

  int halfway[DIM];
  for (int dir = 0;dir < DIM;dir++)
    {
      halfway[dir] = hi[0]/8;
    }
  float eigenvalue = LOfPhi[halfway]/phi[halfway];
  float maxerr=0;
  int imax;
  int tuple[DIM];
  int tuplemax[DIM];
  for (int i = 0;i < domain.sizeOf();i++)
    {
      domain.tupleIndex(i,tuple);
     
      float err =  fabs(LOfPhi[i] - eigenvalue * phi[tuple]);
      if (maxerr < err)
        {
          maxerr =  err;
          imax = i;
          tuplemax[0] = tuple[0];
          tuplemax[1] = tuple[1];
        }
    }
  printf("integer location of error in eigenfunction condition = %d \n", imax);
  printf("tuple locaton of error %d %d \n",tuplemax[0],tuplemax[1]);
  printf("error in eigenfunction condition = %le \n ",maxerr);

  printf("eigenvalue = %le \n ",eigenvalue);
  printf("eigenvalue error = %le \n",eigenvalue - (4*cos(2*M_PI*h) - 4)/h/h);
  //float errexact = 
  printf("exact eigenvalue = %le \n",(4*cos(2*M_PI*h) - 4)/h/h);
}
