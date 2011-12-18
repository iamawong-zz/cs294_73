  
#include "PoissonSolver.H"

PoissonSolver::PoissonSolver()
{
}
PoissonSolver::PoissonSolver(std::tr1::shared_ptr<FFT1D> a_fft1dPtr)
{
  
  m_fft1dptr = a_fft1dPtr;
  m_N = a_fft1dPtr->getN();
  m_M = a_fft1dPtr->getM();
}

void PoissonSolver::solve(MDArray<Real>& a_Rhs) const
{
	//create a contructor
  Box bx = a_Rhs.getBox(); 
  MDArray<complex<double> > fftwForward(bx);
  MDArray<complex<double> > fourierCoef(bx);
  
  double h = 1/m_M;

  for (int k = 0; k < bx.sizeOf();k++)
  {
		fftwForward[k] = complex<double>(a_Rhs[k],0.);
  }
  FFTMD fftmd = FFTMD(m_fft1dptr);
  fftmd.forwardCC(fftwForward);
 
	 for (int i= 0; i<m_M; i++)
	 {
		 for (int j= 0; j<m_M; j++)
		 {
			 int index[2] ={i,j};
			 
			 if(i ==0 && j ==0)
			 {
				 fourierCoef[index] =complex<double>(0,0);
			 }
			 else
			 {
				 complex <double> div(2*cos(2*M_PI*i*h) - 2*cos(2*M_PI*j*h - 4),0);
				 fourierCoef[index] = fftwForward[index]/div;
			 }

		 }
	 }


	 fftmd.inverseCC(fourierCoef);

	 for (int k= 0; k < bx.sizeOf(); k++)
	 {
		
           a_Rhs[k] = real(fourierCoef[k])/m_N;
		 
	 }

}

///  int m_M,m_N; Box m_grid;std::tr1::shared_ptr<FFT1D> m_fft1dptr;
