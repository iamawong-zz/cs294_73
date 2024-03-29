
#ifndef _RK4_H_
#define _RK4_H_

///  generic explicit RK4 algorithm
/**  class X is the type of your state data. class F provides an evaluation
     of the X' for a specified state of X and the time.  If your X' is 
     determined entirely from the current state values then F can ignore.
     F implements the function:
     void operator(X& k, double a_time, double a_dt, const X& a_state, const X& a_shift);

     F interprets the real state as a_state+a_shift

*/
template <class X, class F, class dX>
class RK4
{
public:
  void advance(double a_time, double a_dt, X& a_state);
protected:
  dX m_k; 
  dX m_delta;
  F m_f;
};

//==============================================



template <class X, class F, class dX>
void RK4<X, F, dX>::advance(double a_time, double a_dt, X& a_state)
{
  // this looks odd, but it is a common programming error
  //  you can't write "1/6* m_k"; and get what you want.  the first
  // term is evaluated as integer division, resulting in a zero.
  double sixth = 1, third=1, half = 1;
  sixth/=6; third/=3; half/=2;
  
  m_delta.init(a_state);
  m_k.init(a_state);
  m_f(m_k, a_time, a_dt, a_state, m_k);  // compute k1
  m_delta.increment(sixth, m_k);
  m_k*=half;
  m_f(m_k, a_time+half*a_dt, a_dt, a_state, m_k); // compute k2
  m_delta.increment(third, m_k);
  m_k*=half;
  m_f(m_k, a_time+half*a_dt, a_dt, a_state, m_k);  // conpute k3
  m_delta.increment(third, m_k);
  m_f(m_k, a_time+a_dt, a_dt, a_state, m_k); // compute k4
  m_delta.increment(sixth, m_k);
  a_state.increment(m_delta);
}
  
#endif 
  
