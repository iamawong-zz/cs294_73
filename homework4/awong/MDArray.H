#ifndef _MDARRAY_H_
#define _MDARRAY_H_
#include "Box.H"
template <class T> class MDArray 
{
public:
  /// Default constructor.
  MDArray();
  /// Constructs an MDArray over the Box a_Box;
  MDArray(Box a_box);
  /// Defines a default-constructed MDArray.  
  void define(Box a_box);
  /// Destructor.
  ~MDArray(); 
  /// Setting an MDArray to a constant value.
  void setVal(T a_val);
  /// Copying the contents to the input MDArray on the intersection of the boxes.
  void copyTo(MDArray<T>& a_destArray) const;
  /// Pointwise addition operator. a_rhs.m_box must be equal to this->m_box.
  void operator+=(const MDArray& a_rhs);
  /// pointwise subtraction operator. a_rhs.m_box must be equal to this->m_box.
  void operator-=(const MDArray& a_rhs); 
  /// Pointwise multiplication operator. a_rhs.m_box must be equal to this->m_box.
  void operator*=(const MDArray& a_rhs);
  /// Pointwise division operator. a_rhs.m_box must be equal to this->m_box.
  void operator/=(const MDArray& a_rhs);
  
  /// Shift operators.
  /// Should be void.
  void shift(int a_tuple[DIM]);
  void shift(int a_dir, int a_len);
  /// Indexing operator.
  T& operator[](int a_index[DIM]) const;
  /// Indexing shifted operator.
  T& indexShift(const int a_tuple[DIM], const int& a_dir,const int& a_shift) const;
  /// Indexing into 1D array representation of data.
  T& operator[](int a_linearIndex) const;
  /// Access the box.
  const Box& getBox() const;
  ///Add a T to all elements.
  void operator+=(const T& a_rhs);
  /// Subtract a T from all elements.
  void operator-=(const T& a_rhs);
  /// Multiply all elements by a scalar.
  void operator*=(const T& a_rhs);
  /// Divide all elements by a scalar.
  void operator/=(const T& a_rhs);
private:
  T *m_data;
  Box m_box;
  bool m_isDefined;
  ///we make copy constructors, assignment operators private to prevent their use.
  MDArray(const MDArray& a_MDArray);
  MDArray& operator=(const MDArray& a_MDArray);
};

#include "MDArrayImplem.H"  
#endif
