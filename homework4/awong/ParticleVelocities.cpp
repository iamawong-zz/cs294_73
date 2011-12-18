#include "ParticleVelocities.H"

void ParticleVelocities::operator()(ParticleShift& a_kOut, 
		double a_time, double dt, 
		const ParticleSet& a_state,
		const ParticleShift& a_kIn)
{
  int binIndex[DIM];
  int highCorner[DIM];
  int binTuple[DIM];
  ParticleSet Splus = a_state;
  /// Incrementing with K
  Splus.increment(a_kIn);
  int N = Splus.m_box.sizeOf();
  MLCVortexMethod vortex(Splus.m_cutoff, Splus.m_dx, Splus.m_lowCorner);
  MDArray<Real> totalVelocities[DIM];
  MDArray<double> RHSlargex(Splus.m_box);  
  MDArray<double> RHSlargey(Splus.m_box);

  Splus.m_box.getHighCorner(highCorner);

  /// This is to find how many intervals on each side
  for(unsigned int i=0; i<DIM; i++) {
    binIndex[i] = Splus.m_box.length(i)/Splus.m_dx;
    totalVelocities[i].define(Splus.m_box);
    totalVelocities[i].setVal(0);
  }

  PoissonInfinite pSolver(binIndex[0], 1/binIndex[0], binIndex[0]/8);

  /// Rebinning
  Splus.rebin(Splus.m_box, Splus.m_dx, Splus.m_lowCorner);
  
  /// Start going through the bins
  for(int i=0; i<N; i++) {
    if(!(*Splus.m_bins)[i].empty()) {
      MDArray<Real> gridVelocities[DIM];
      Splus.m_box.tupleIndex(i, binTuple);
     
      int low[DIM];
      int high[DIM];
      low[0] = binTuple[0]-Splus.m_cutoff-1;
      low[1] = binTuple[1]-Splus.m_cutoff-1;
      high[0] = binTuple[0]+Splus.m_cutoff+1;
      high[1] = binTuple[1]+Splus.m_cutoff+1; 
	
      /// This is the box of C+1
      Box localBox(low, high);
      for(int i=0; i<DIM; i++) {
	gridVelocities[i].define(localBox);
	gridVelocities[i].setVal(0);
      }
      
      /// Boxes of C
      MDArray<double> RHSx(localBox.grow(-1));  
      MDArray<double> RHSy(localBox.grow(-1));
      
      vortex.computeVelocities(gridVelocities, Splus, binTuple);
      pSolver.applyOperator(RHSx, gridVelocities[0], localBox);
      pSolver.applyOperator(RHSy, gridVelocities[1], localBox); 
      RHSlargex += RHSx;
      RHSlargey += RHSy;
    }
  }/// Finish going through the bins
  /// Solving over the whole grid
  pSolver.solve(totalVelocities[0], RHSlargex);
  pSolver.solve(totalVelocities[1], RHSlargey);
  for(int i=0; i<N; i++) {
    if(!(*Splus.m_bins)[i].empty()) {
      Splus.m_box.tupleIndex(i, binTuple);
      vortex.MLCInterpolate(a_kOut, Splus, totalVelocities, binTuple);
    }
  }/// Finish going through the bins  
  a_kOut *= dt;
}
