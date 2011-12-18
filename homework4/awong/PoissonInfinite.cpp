#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "MDArray.H"
#include "FFTW1DSine.H"
#include "FFTMDSine.H"
#include "PoissonInfinite.H"
using namespace std;
PoissonInfinite::PoissonInfinite()
{
  assert(DIM==2);
};
PoissonInfinite::PoissonInfinite(unsigned int a_N, double a_h, unsigned int a_buffersize)
{
  assert(DIM == 2);
  m_N = a_N; 
  m_L = m_N*a_h;
  m_fftmd = new FFTMDSine(a_N-1);
  m_fftmdBig = new FFTMDSine(a_N + 2*a_buffersize-1);
  m_h = a_h;
  int lo[DIM]={0,0};
  int hi[DIM]={m_N,m_N};
  m_box = Box(lo,hi);
  m_boxBig = m_box.grow(a_buffersize);
  int weights[9] = { 20,-4,-4,-4,-4,-1,-1,-1,-1};
  for (int l  = 0; l < STENCILSIZE;l++)
    {
      m_weights[l] = weights[l];
    }
  //int offsets[STENCILSIZE][DIM] ={0,0,1,0,0,1,-1,0,0,-1,1,1,-1,1,1,-1,1,-1};
  int offsets[STENCILSIZE][DIM] ={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{1,-1},{-1,-1}}; 
   for (int l = 0;l < STENCILSIZE;l++)
     {
       m_offsets[l][0] = offsets[l][0];
       m_offsets[l][1] = offsets[l][1];
       } 
};
void
PoissonInfinite::solve(MDArray<double>& a_solution, const MDArray<double>& a_rhs) const
{
  MDArray<double> bigBoundaryCharges[DIM][2];
  MDArray<double> boundaryCharges[DIM][2];
  MDArray<double> copy(m_box);
  a_rhs.copyTo(copy);
  dirichletSolve(copy,true);
  computeInnerBoundaryCharges(boundaryCharges,copy);
  MDArray<double> bigField(m_boxBig);
  computeOuterBoundaryCharges(bigBoundaryCharges,boundaryCharges);
 
  bigField.setVal(0.);
  a_rhs.copyTo(bigField); 
  /* for (int dir = 0; dir < DIM;dir++)
    { 
      for (int lohi = 0; lohi < 2;lohi++)
        {
          bigBoundaryCharges[dir][lohi].copyTo(bigField);
        }
        }*/
  dirichletSolve(bigField,false);
  bigField.copyTo(a_solution);
};
void 
PoissonInfinite::dirichletSolve(MDArray<double>& a_solution,bool a_isInner) const
{
  int lowCorner[DIM];
  int highCorner[DIM];
  int NptsLogical;
  a_solution.getBox().getLowCorner(lowCorner);  
  a_solution.getBox().getHighCorner(highCorner);
  NptsLogical = highCorner[0] - lowCorner[0];
  for (int dir = 0;dir < DIM; dir++)
    {
      lowCorner[dir] *=-1;
    }
  Box bx;
  a_solution.shift(lowCorner);
  bx = a_solution.getBox();
  if (a_isInner) 
    {
      m_fftmd->transform(a_solution);
    }
  else
    {
      m_fftmdBig->transform(a_solution);
    }
  int lowTrans[DIM] = {1,1};
  int highTrans[DIM];
  bx.getHighCorner(highTrans);
  Box bxTrans = Box(lowTrans,highTrans);

   for (int k=0;k < bxTrans.sizeOf();k++)
    {
      int tuple[DIM];
      
      bxTrans.tupleIndex(k,tuple);
      
      double symbol = 
        (20. - 
         8*cos(M_PI*tuple[0]/NptsLogical) -
         8*cos(M_PI*tuple[1]/NptsLogical) -
         4*cos(M_PI*tuple[0]/NptsLogical)*cos(M_PI*tuple[1]/NptsLogical))
        /m_h/m_h/6;
      a_solution[tuple] = a_solution[tuple]/symbol;
      }
  if (a_isInner) 
    {
      m_fftmd->transform(a_solution);
    }
  else
    {
      m_fftmdBig->transform(a_solution);
    }
  for (int dir = 0;dir < DIM; dir++)
    {
      lowCorner[dir] *=-1;
    }
  a_solution.shift(lowCorner);
}
void 
PoissonInfinite::computeInnerBoundaryCharges(MDArray<double> a_charges[DIM][2],
                                             const MDArray<double>& a_field) const
{
  int lowCorner[DIM],highCorner[DIM];  
  a_field.getBox().getLowCorner(lowCorner);
  a_field.getBox().getHighCorner(highCorner);
  int lowBoundary[DIM],highBoundary[DIM];
  int tuple[DIM];
  int tuple0[DIM];
  int tuplep[DIM];
  int tuplem[DIM];
  // Make containers.
  Box bx[DIM][2];
  
  for (int dir = 0; dir <DIM ; dir++)
    {
      // Low boundary charge container.
      lowBoundary[dir] = lowCorner [dir];
      highBoundary[dir] = lowBoundary[dir];
      lowBoundary[(dir+1)%DIM] = lowCorner[(dir+1)%DIM]  + dir;
      highBoundary[(dir+1)%DIM] = highCorner[(dir+1)%DIM] - dir;
      
      Box bb = Box(lowBoundary,highBoundary);
      a_charges[dir][0].define(bb);
      bx[dir][0] = bb;
      

      // high boundary charge container/
      
      highBoundary[dir] = highCorner[dir];
      lowBoundary[dir] = highCorner[dir];
      bb = Box(lowBoundary,highBoundary);
      a_charges[dir][1].define(bb);
      bx[dir][1] = bb;
    }
  
  for (int dir = 0; dir < DIM; dir++)
    {
      for (int lohi = 0; lohi < 2; lohi++)
        {
          MDArray<double>& charge = a_charges[dir][lohi];
          Box bbField = bx[dir][lohi].grow(1);
          MDArray<double> field(bbField);
          field.setVal(0.);
          a_field.copyTo(field);
          applyOperator(charge,field,bx[dir][lohi]);
        }
    }   
};
void 
PoissonInfinite::computeOuterBoundaryCharges(
                     MDArray<double> a_outerCharges[DIM][2],
               const MDArray<double> a_innerCharges[DIM][2]) const
{
int lowCorner[DIM],highCorner[DIM];  
  m_boxBig.getLowCorner(lowCorner);
  m_boxBig.getHighCorner(highCorner);
  int lowBoundary[DIM],highBoundary[DIM];
  int tuple[DIM]; 
  Box outerBoundary[DIM][2];
  
  // Make containers.
  Box bx[DIM][2];

  for (int dir = 0; dir <DIM ; dir++)
    {
      // Low boundary charge container.
      lowBoundary[dir] = lowCorner[dir]+1;
      highBoundary[dir] = lowBoundary[dir];
      lowBoundary[(dir+1)%DIM] = lowCorner[(dir+1)%DIM] + 1;
      highBoundary[(dir+1)%DIM] = highCorner[(dir+1)%DIM] - 1;
      
      Box bb = Box(lowBoundary,highBoundary);
      a_outerCharges[dir][0].define(bb);
      a_outerCharges[dir][0].setVal(0.);
      bb = bb.shift(dir,-1);
      bb = bb.grow((dir + 1)%DIM,1-dir);
      outerBoundary[dir][0] = bb;

      // high boundary charge container/
      
      highBoundary[dir] = highCorner[dir] - 1;
      lowBoundary[dir] = highCorner[dir] - 1;
      bb = Box(lowBoundary,highBoundary);
      a_outerCharges[dir][1].define(bb);
      a_outerCharges[dir][1].setVal(0.);   
      bb = bb.shift(dir,1);
      bb = bb.grow((dir + 1)%DIM,1-dir);
      outerBoundary[dir][1] = bb;
     
    }
  for (int dir =0;dir < DIM;dir++)
    {
      for (int hiLo = 0;hiLo < 2; hiLo++)
        {
          for (int dir2 = 0;dir2 < DIM; dir2++)
            {
              for (int hiLo2 = 0; hiLo2 < 2;hiLo2++)
                {
                  faceToFace(a_outerCharges[dir][hiLo],
                             a_innerCharges[dir2][hiLo2],
                             outerBoundary[dir][hiLo]);
                }
            }
        }
    }
};
void
PoissonInfinite::faceToFace(MDArray<double>& a_outerCharge,
                            const MDArray<double>& a_innerCharge,
                            const Box& a_outerBoundary) const
{
  int tuple[DIM];
  int tuple2[DIM];
  int lowCorner[DIM];
  int highCorner[DIM];
  int lowBoundary[DIM];
  int highBoundary[DIM];
 
  Box fieldBox = a_outerCharge.getBox().grow(1);
  Box innerChargeBox = a_innerCharge.getBox();
  MDArray<double> outerField(fieldBox);
  outerField.setVal(0.);
  for (int k = 0;k < innerChargeBox.sizeOf();k++)
    {
      innerChargeBox.tupleIndex(k,tuple);
      for (int k2 = 0; k2 < a_outerBoundary.sizeOf();k2++)
        {
          a_outerBoundary.tupleIndex(k2,tuple2);
          double mag = 0;
          for (int dir = 0;dir < DIM;dir++)
            {
              int dtuple = tuple2[dir] - tuple[dir];
              mag += pow(dtuple*m_h,2);
            }
          outerField[tuple2] += -m_h*m_h*a_innerCharge[tuple]*log(mag)/4/M_PI;
        }
    }
  Box bxCharge = a_outerCharge.getBox();
  MDArray<double> charge(bxCharge);
  applyOperator(charge,outerField,bxCharge);
  a_outerCharge += charge; 
}
void PoissonInfinite::applyOperator(MDArray<double>& a_LOfPhi, 
                                    const MDArray<double>& a_phi,
                                    const Box& a_bxApply) const
{
    int stencil[STENCILSIZE][DIM];
    int tuple[DIM];
    
    for (int k = 0;k < a_bxApply.sizeOf();k++)
      {
        a_bxApply.tupleIndex(k,tuple);
        getStencil(stencil,tuple);
        a_LOfPhi[tuple] = 0.;
        for (int l =0;l < STENCILSIZE;l++)
          {
            a_LOfPhi[tuple] += a_phi[stencil[l]]*m_weights[l];
          }
        a_LOfPhi[tuple]/= 6*m_h*m_h;
      }
   
}
void PoissonInfinite::getStencil(int a_stencil[STENCILSIZE][DIM],int tuple[DIM]) const
{
  for (int l = 0;l < STENCILSIZE;l++)
    {
      for (int dir = 0;dir < DIM; dir++)
        {
          a_stencil[l][dir] = tuple[dir] + m_offsets[l][dir];
        }
    }
}
