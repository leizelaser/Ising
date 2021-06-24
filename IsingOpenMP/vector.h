/*
  Basic vector library
*/
#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <fstream>
#include <math.h>

#define DIM 3

const unsigned PDIM = unsigned(pow(3,DIM));
const unsigned DDIM = unsigned(pow(2,DIM));
const unsigned PIDDIM = unsigned(pow(2,2*DIM));

template <unsigned D, typename T=double>
class vector {
  
  public:
  T x[D];
  
  public:
  
  vector();
  vector(const T);
  vector(const T[D]);           
  vector(const vector&);
  ~vector();
  void setzero();

  //Scalar operations
  vector<D, T>& operator+=(const T);
  vector<D, T>& operator-=(const T);
  vector<D, T>& operator*=(const T);
  vector<D, T>& operator/=(const T);
  vector<D, T>& operator%=(const T);


  vector<D, T> operator+(const T) const;
  vector<D, T> operator-(const T) const;
  vector<D, T> operator*(const T) const;
  vector<D, T> operator/(const T) const;
  vector<D, T> operator%(const T) const;

  //Vectorized operations
  
  vector<D, T>& operator+=(const vector<D, T>&);
  vector<D, T>& operator-=(const vector<D, T>&);
  vector<D, T>& operator*=(const vector<D, T>&);
  vector<D, T>& operator/=(const vector<D, T>&);
  vector<D, T>& operator%=(const vector<D, T>&);

  vector<D, T> operator+(const vector<D, T>&) const;
  vector<D, T> operator-(const vector<D, T>&) const;
  vector<D, T> operator*(const vector<D, T>&) const;
  vector<D, T> operator/(const vector<D, T>&) const;
  vector<D, T> operator%(const vector<D, T>&) const;

  //Simple Self-Modulo (plus or minus 1)
  vector<D, T>& selfmodp(const T);
  vector<D, T>& selfmodp(const vector<D, T>&);

  vector<D, T>& selfmodm(const T);
  vector<D, T>& selfmodm(const vector<D, T>&);

  //Boolean operator
  bool operator==(const vector<D, T> &a) const;

  //Vectorized Boolean
  vector<D, bool> operator<=(const vector<D, T> &a) const;
  vector<D, bool> operator>=(const vector<D, T> &a) const;
  vector<D, bool> operator<(const vector<D, T> &a) const;
  vector<D, bool> operator>(const vector<D, T> &a) const;

  vector<D, bool> operator<=(const double a) const;
  vector<D, bool> operator>=(const double a) const;
  vector<D, bool> operator<(const double a) const;
  vector<D, bool> operator>(const double a) const;

  //Type transform
  vector<D, int> Integer() const;
  vector<D, unsigned> Unsigned() const;
  vector<D, double> Double() const;

  static vector<D, int> Integer(const vector<D, T>&); 
  static vector<D, unsigned> Unsigned(const vector<D, T>&);
  static vector<D, double> Double(const vector<D, T>&); 

  //Ceiling, Floor
  vector<D, int> Ceiling() const;
  static vector<D, int> Ceiling(const vector<D, T>&); 

  vector<D, int> FloorInt() const;
  static vector<D, int> FloorInt(const vector<D, T>&); 
  
  vector<D, double> FloorDouble() const;
  static vector<D, double> FloorDouble(const vector<D, T>&); 
  
  //Access
  T& operator[](const unsigned);
 
  //Inner product
  double dot(const vector<D, T>&) const;
  static double dot(const vector<D, T>&, const vector<D, T>&);
  
  double norm_squared() const;
  static double norm_squared(const vector<D, T>&);
  double norm() const;
  
  //Statistics
  T maximum() const;
  static T maximum(const vector<D, T>&);
  T minimum() const;
  static T minimum(const vector<D, T>&);

  bool AnyT(T) const;
  static bool AnyT(const vector<D, T>&, T);
  bool AllT(T) const;
  static bool AllT(const vector<D, T>&, T);

  bool AnyV(const vector<D, T>&) const;
  static bool AnyV(const vector<D, T>&, const vector<D, T>&);
  bool AllV(const vector<D, T>&) const;
  static bool AllV(const vector<D, T>&, const vector<D, T>&);

  T sum() const;
  static T sum(const vector<D, T>&);
  
  //Shuffle
  void shuffle(vector<D, T>&);
  vector<D, double>& SetRandDouble(); 
  static vector<D, double>& SetRandDouble(vector<D, double>&); 
  
  //Read & write
  void read(std::ifstream&);
  void write(std::ofstream&) const;
};

template <unsigned D, typename T>
std::ostream& operator<<(std::ostream&, const vector<D, T>&);

// constructor
// ~~~~~~~~~~~
template <unsigned D, typename T>
vector<D, T>::vector()
{
  for(unsigned k=0; k<D; k++)
    x[k] = 0;
}


template <unsigned D, typename T>
vector<D, T>::vector(const T x_i)
{
  for(unsigned k=0; k<D; k++)
    x[k] = x_i;
}


template <unsigned D, typename T>
vector<D, T>::vector(const T x_i[D])
{
  for(unsigned k=0; k<D; k++)
    x[k] = x_i[k];
}


template <unsigned D, typename T>
vector<D, T>::vector(const vector<D, T> &v)
{
  for(unsigned k=0; k<D; k++)
    x[k] = v.x[k];
}


// destructor
// ~~~~~~~~~~
template <unsigned D, typename T>
vector<D, T>::~vector()
{
}

// +=(scalar)
// ~~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T>& vector<D, T>::operator+=(const T s)
{
  for(unsigned k=0; k<D; k++)
    x[k] += s;
  return *this;
}

// +=(vector)
// ~~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T>& vector<D, T>::operator+=(const vector<D, T> &v)
{
  for(unsigned k=0; k<D; k++)
    x[k] += v.x[k];
  return *this;
}

// -=(scalar)
// ~~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T>& vector<D, T>::operator-=(const T s)
{
  for(unsigned k=0; k<D; k++)
    x[k] -= s;
  return *this;
}

// -=(vector)
// ~~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T>& vector<D, T>::operator-=(const vector<D, T> &v)
{
  for(unsigned k=0; k<D; k++)
    x[k] -= v.x[k];
  return *this;
}

// *=(scalar)
// ~~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T>& vector<D, T>::operator*=(const T s)
{
  for(unsigned k=0; k<D; k++)
    x[k] *= s;
  return *this;
}

// *=(vector)
// ~~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T>& vector<D, T>::operator*=(const vector<D, T> &v)
{
  for(unsigned k=0; k<D; k++)
    x[k] *= v.x[k];
  return *this;
}

// /=(scalar)
// ~~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T>& vector<D, T>::operator/=(const T s)
{
  for(unsigned k=0; k<D; k++)
    x[k] /= s;
  return *this;
}

// /=(vector)
// ~~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T>& vector<D, T>::operator/=(const vector<D, T> &v)
{
  for(unsigned k=0; k<D; k++)
    x[k] /= v.x[k];
  return *this;
}

// %=(scalar, postive)
// ~~~~~~~~~~~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T>& vector<D, T>::operator%=(const T s)
{
  for(unsigned k=0; k<D; k++)
  {
    x[k] = x[k] % s;
    if(x[k] < 0)
      x[k] += s;
  }
  return *this;
}

// %=(vector, postive)
// ~~~~~~~~~~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T>& vector<D, T>::operator%=(const vector<D, T> &a)
{
  for(unsigned k=0; k<D; k++)
  {
    x[k] = x[k] % a.x[k];
    if(x[k] < 0)
      x[k] += a.x[k];
  }
  return *this;
}


// +(scalar)
// ~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T> vector<D, T>::operator+(const T s) const
{
  vector<D, T> c;
  for(unsigned k=0; k<D; k++)
    c.x[k] = x[k] + s;
  return c;
}

// +(vector)
// ~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T> vector<D, T>::operator+(const vector<D, T> &a) const
{
  vector<D, T> c;
  for(unsigned k=0; k<D; k++)
    c.x[k] = x[k] + a.x[k];
  return c;
}


// -(scalar)
// ~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T> vector<D, T>::operator-(const T s) const
{
  vector<D, T> c;
  for(unsigned k=0; k<D; k++)
    c.x[k] = x[k] - s;
  return c;
}

// -(vector)
// ~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T> vector<D, T>::operator-(const vector<D, T> &a) const
{
  vector<D, T> c;
  for(unsigned k=0; k<D; k++)
    c.x[k] = x[k] - a.x[k];
  return c;
}

// *(scalar)
// ~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T> vector<D, T>::operator*(const T s) const
{
  vector<D, T> c;
  for(unsigned k=0; k<D; k++)
    c.x[k] = x[k] * s;
  return c;
}

// *(vector)
// ~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T> vector<D, T>::operator*(const vector<D, T> &a) const
{
  vector<D, T> c;
  for(unsigned k=0; k<D; k++)
    c.x[k] = x[k] * a.x[k];
  return c;
}

// /(scalar)
// ~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T> vector<D, T>::operator/(const T s) const
{
  vector<D, T> c;
  for(unsigned k=0; k<D; k++)
    c.x[k] = x[k] / s;
  return c;
}

// /(vector)
// ~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T> vector<D, T>::operator/(const vector<D, T> &a) const
{
  vector<D, T> c;
  for(unsigned k=0; k<D; k++)
    c.x[k] = x[k] / a.x[k];
  return c;
}

// %(scalar, postive)
// ~~~~~~~~~~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T> vector<D, T>::operator%(const T s) const
{
  vector<D, T> c;
  for(unsigned k=0; k<D; k++)
  {
    c.x[k] = x[k] % s;
    if(c.x[k] < 0)
      c.x[k] += s;
  }
  return c;
}

// %(vector, postive)
// ~~~~~~~~~~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T> vector<D, T>::operator%(const vector<D, T> &a) const
{
  vector<D, T> c;
  for(unsigned k=0; k<D; k++)
  {
    c.x[k] = x[k] % a.x[k];
    if(c.x[k] < 0)
      c.x[k] += a.x[k];
  }
  return c;
}

// Simple Modulo (plus 1 if negative)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T>& vector<D, T>::selfmodp(const T s)
{
  for(unsigned k=0; k<D; k++)
  {
    if(x[k] < 0)
      x[k] += s;
  }
  return *this;
}

template <unsigned D, typename T>
inline vector<D, T>& vector<D, T>::selfmodp(const vector<D, T> &a)
{
  for(unsigned k=0; k<D; k++)
  {
    if(x[k] < 0)
      x[k] += a.x[k];
  }
  return *this;
}

// Simple Modulo (minus 1 if >1)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <unsigned D, typename T>
inline vector<D, T>& vector<D, T>::selfmodm(const T s)
{
  for(unsigned k=0; k<D; k++)
  {
    if(x[k] >= s)
      x[k] -= s;
  }
  return *this;
}

template <unsigned D, typename T>
inline vector<D, T>& vector<D, T>::selfmodm(const vector<D, T> &a)
{
  for(unsigned k=0; k<D; k++)
  {
    if(x[k] >= a.x[k])
      x[k] -= a.x[k];
  }
  return *this;
}


// ==
// ~~
template <unsigned D, typename T>
inline bool vector<D, T>::operator==(const vector<D, T> &a) const
{
  for(unsigned k=0; k<D; k++)
    {
      if (!(x[k]==a.x[k]))
	return false;
    }
  return true;
}

// <=
// ~~
template <unsigned D, typename T>
inline vector<D, bool> vector<D, T>::operator<=(const vector<D, T> &a) const
{
  vector<D, bool> b;
  for(unsigned k=0; k<D; k++)
    b[k] = (x[k] <= a.x[k]);
  return b;
}

template <unsigned D, typename T>
inline vector<D, bool> vector<D, T>::operator<=(const double a) const
{
  vector<D, bool> b;
  for(unsigned k=0; k<D; k++)
    b[k] = (x[k] <= a);
  return b;
}

// >=
// ~~
template <unsigned D, typename T>
inline vector<D, bool> vector<D, T>::operator>=(const vector<D, T> &a) const
{
  vector<D, bool> b;
  for(unsigned k=0; k<D; k++)
    b[k] = (x[k] >= a.x[k]);
  return b;
}

template <unsigned D, typename T>
inline vector<D, bool> vector<D, T>::operator>=(const double a) const
{
  vector<D, bool> b;
  for(unsigned k=0; k<D; k++)
    b[k] = (x[k] >= a);
  return b;
}

// <
// ~
template <unsigned D, typename T>
inline vector<D, bool> vector<D, T>::operator<(const vector<D, T> &a) const
{
  vector<D, bool> b;
  for(unsigned k=0; k<D; k++)
    b[k] = (x[k] < a.x[k]);
  return b;
}

template <unsigned D, typename T>
inline vector<D, bool> vector<D, T>::operator<(const double a) const
{
  vector<D, bool> b;
  for(unsigned k=0; k<D; k++)
    b[k] = (x[k] < a);
  return b;
}
// >
// ~
template <unsigned D, typename T>
inline vector<D, bool> vector<D, T>::operator>(const vector<D, T> &a) const
{
  vector<D, bool> b;
  for(unsigned k=0; k<D; k++)
    b[k] = (x[k] > a.x[k]);
  return b;
}

template <unsigned D, typename T>
inline vector<D, bool> vector<D, T>::operator>(const double a) const
{
  vector<D, bool> b;
  for(unsigned k=0; k<D; k++)
    b[k] = (x[k] > a);
  return b;
}


// Integer
// ~~~~~~~
template <unsigned D, typename T>
inline vector<D, int> vector<D, T>::Integer() const
{
  vector<D, int> c;

  for(unsigned k=0; k<D; k++)
    c[k] = (int)x[k];
  return c;
}

template <unsigned D, typename T>
inline vector<D, int> vector<D, T>::Integer(const vector<D, T>& v)
{
  return v.Integer();
}

// Unsigned
// ~~~~~~~~

template <unsigned D, typename T>
inline vector<D, unsigned> vector<D, T>::Unsigned() const
{
  vector<D, unsigned> c;

  for(unsigned k=0; k<D; k++)
    c[k] = (unsigned)x[k];
  return c;
}

template <unsigned D, typename T>
inline vector<D, unsigned> vector<D, T>::Unsigned(const vector<D, T>& v)
{
  return v.Unsigned();
}

// Double
// ~~~~~~~
template <unsigned D, typename T>
inline vector<D, double> vector<D, T>::Double() const
{
  vector<D, double> c;

  for(unsigned k=0; k<D; k++)
    c[k] = (double)x[k];
  return c;
}

template <unsigned D, typename T>
inline vector<D, double> vector<D, T>::Double(const vector<D, T>& v)
{
  return v.Double();
}

// Ceiling
template <unsigned D, typename T>
inline vector<D, int> vector<D, T>::Ceiling() const
{
  vector<D, int> c;
  for(unsigned k=0; k<D; k++)
  {
    c[k] = (int)x[k];
    c[k] += (x[k] - c[k] > 0 ? 1:0);
  }  
  return c;
}

template <unsigned D, typename T>
inline vector<D, int> vector<D, T>::Ceiling(const vector<D, T>& v) 
{
  return v.Ceiling();
}

// Floor
template <unsigned D, typename T>
inline vector<D, int> vector<D, T>::FloorInt() const
{
  vector<D, int> c;
  for(unsigned k=0; k<D; k++)
  {
    c[k] = (int)x[k];
    c[k] += (x[k] - c[k] < 0 ? (-1):0);
  }  
  return c;
}

template <unsigned D, typename T>
inline vector<D, int> vector<D, T>::FloorInt(const vector<D, T>& v) 
{
  return v.FloorInt();
}

template <unsigned D, typename T>
inline vector<D, double> vector<D, T>::FloorDouble() const
{
  vector<D, double> c;
  for(unsigned k=0; k<D; k++)
  {
    c[k] = double((int)x[k]);
    c[k] += (x[k] - c[k] < 0 ? (-1.0):0.0);
  }  
  return c;
}

template <unsigned D, typename T>
inline vector<D, double> vector<D, T>::FloorDouble(const vector<D, T>& v) 
{
  return v.FloorDouble();
}

// []
// ~~
template <unsigned D, typename T>
inline T& vector<D, T>::operator[](const unsigned i)
{
  return x[i];
}


// Dot
// ~~~
template <unsigned D, typename T>
inline double vector<D, T>::dot(const vector<D, T> &a) const
{
  double d=0;

  for(unsigned k=0; k<D; k++)
    d += x[k] * a.x[k];

  return d;
}

template <unsigned D, typename T>
inline double vector<D, T>::dot(const vector<D, T> &a, const vector<D, T> &b)
{
  return a.dot(b);
}


// NormSquared
// ~~~~~~~~~~~
template <unsigned D, typename T>
inline double vector<D, T>::norm_squared() const
{
  return dot(*this, *this);
}

// Norm
// ~~~~~~~~~~~
template <unsigned D, typename T>
inline double vector<D, T>::norm() const
{
  return sqrt(dot(*this, *this));
}


template <unsigned D, typename T>
inline double vector<D, T>::norm_squared(const vector<D, T>& v)
{
  return v.norm_squared();
}

// Maximum & minimum
// ~~~~~~~~~~~~~~~~~
template <unsigned D, typename T>
inline T vector<D, T>::maximum() const
{
  return maximum(*this);
}

template <unsigned D, typename T>
inline T vector<D, T>::maximum(const vector<D, T> &a)
{
  T t = a.x[0];
  for (unsigned i=1; i<D; i++)
  {
    if (t > a.x[i])
      t = a.x[i];
  }
  return t;
}

template <unsigned D, typename T>
inline T vector<D, T>::minimum() const
{
  return minimum(*this);
}

template <unsigned D, typename T>
inline T vector<D, T>::minimum(const vector<D, T> &a)
{
  T t = a.x[0];
  for (unsigned i=1; i<D; i++)
  {
    if (t < a.x[i])
      t = a.x[i];
  }
  return t;
}

// Any & All
// ~~~~~~~~~

template <unsigned D, typename T>
inline bool vector<D, T>::AnyT(T t) const
{
  return AnyT(*this, t);
}

template <unsigned D, typename T>
inline bool vector<D, T>::AnyT(const vector<D, T> &a, T t)
{
  for (unsigned i=0; i<D; i++)
  {
    if (t == a.x[i])
      return true;
  }
  return false;
}

template <unsigned D, typename T>
inline bool vector<D, T>::AllT(T t) const
{
  return AllT(*this, t);
}

template <unsigned D, typename T>
inline bool vector<D, T>::AllT(const vector<D, T> &a, T t)
{
  for (unsigned i=0; i<D; i++)
  {
    if (t != a.x[i])
      return false;
  }
  return true;
}


template <unsigned D, typename T>
inline bool vector<D, T>::AnyV(const vector<D, T> &v) const
{
  return AnyV(*this, v);
}

template <unsigned D, typename T>
inline bool vector<D, T>::AnyV(const vector<D, T> &a, const vector<D, T> &v)
{
  for (unsigned i=0; i<D; i++)
  {
    if (v.x[i] == a.x[i])
      return true;
  }
  return false;
}

template <unsigned D, typename T>
inline bool vector<D, T>::AllV(const vector<D, T> &v) const
{
  return AllV(*this, v);
}

template <unsigned D, typename T>
inline bool vector<D, T>::AllV(const vector<D, T> &a, const vector<D, T> &v)
{
  for (unsigned i=0; i<D; i++)
  {
    if (v.x[i] != a.x[i])
      return false;
  }
  return true;
}

template <unsigned D, typename T>
inline T vector<D, T>::sum() const
{
  return sum(*this);
};
 
template <unsigned D, typename T>
inline T vector<D, T>::sum(const vector<D, T>&a)
{
  T t = a.x[0];
  for (unsigned i=1; i<D; i++)
    t += a.x[i];
  return t;
};

// shuffle
// ~~~~~~~
template <unsigned D, typename T>
inline void vector<D, T>::shuffle(vector<D, T>&)
{
  unsigned d(D-1), swapind;
  T t;
  while (d>1)
  {
    swapind = rand()%D;
    if (swapind == d)
      continue;
    else
    {
      t = x[d];
      x[d] = x[swapind]; 
      x[swapind] = t;
    }
    d--;
  }
}

// set rand
// ~~~~~~~~
template <unsigned D, typename T>
inline vector<D, double>& vector<D, T>::SetRandDouble()
{
  return SetRandDouble(*this);
}

template <unsigned D, typename T>
inline vector<D, double>& vector<D, T>::SetRandDouble(vector<D, double>& v)
{
  for (int i=0; i<D; i++)
    v[i] = double(rand())/RAND_MAX;
  return v;
}

// read
// ~~~~
template <unsigned D, typename T>
void vector<D, T>::read(std::ifstream& in)
{
  in.read((char*)x, sizeof(T)*D);
}

// write
// ~~~~~
template <unsigned D, typename T>
void vector<D, T>::write(std::ofstream& out) const
{
  out.write((const char*)x, sizeof(T)*D);
}

// set all element to zero
// ~~~~~~~~~~~~~~~~~~~~~~~
template <unsigned D, typename T>
inline void vector<D, T>::setzero()
{
  for(unsigned k=0; k<D; k++)
    x[k] = 0;
}

// Insertion
// ~~~~~~~~~
template <unsigned D, typename T>
std::ostream& operator<<(std::ostream& os, const vector<D, T>& v)
{
  os << "(";

  for(unsigned k=0; k<D-1; k++)
    os << v.x[k] << ", ";

  os << v.x[D-1] << ")";

  return os;
}


// ======================================================================
// vector_field 
// ======================================================================

// A field of V-vectors on a D dimensional manifold

template<unsigned V, unsigned D, typename T=double>
class vector_field {
  
 public:
  unsigned elements;

 private:
  vector<V, T>* f;
  vector<D, unsigned> size;           // number of grid points for each dimension
  vector<D, unsigned> offset;
 
 public:

  vector_field();
  vector_field(const vector<D, unsigned>&);
  ~vector_field();

  vector<D, unsigned> get_size() const;
  void set_size(const vector<D, unsigned>&);

  vector<V, T>& get(const vector<D, unsigned>&);

  void read(std::ifstream&);
  void write(std::ofstream&) const;

  static void swap(vector_field<V, D, T>&, vector_field<V, D, T>&);
};


// vector_field
// ~~~~~~~~~~~~
template<unsigned V, unsigned D, typename T>
vector_field<V, D, T>::vector_field()
  : f(0), elements(0)
{
}


// vector_field
// ~~~~~~~~~~~~
template<unsigned V, unsigned D, typename T>
vector_field<V, D, T>::vector_field(const vector<D, unsigned>& s)
  : f(0)
{
  set_size(s);
}

// ~vector_field
// ~~~~~~~~~~~~~
template <unsigned V, unsigned D, typename T>
vector_field<V, D, T>::~vector_field()
{
  if(f != 0)
    delete[] f;
}

// get_size
// ~~~~~~~~
template<unsigned V, unsigned D, typename T>
inline vector<D, unsigned> vector_field<V, D, T>::get_size() const
{
  return size;
}

// set_size
// ~~~~~~~~
template<unsigned V, unsigned D, typename T>
void vector_field<V, D, T>::set_size(const vector<D, unsigned>& s)
{
  if(f != 0)
    delete[] f;

  size = s;

  elements = 1;
  for(unsigned i=0; i<D; i++) {
    offset[i] = elements;
    elements *= size.x[i];
  }

  f = new vector<V, T>[elements];
}

// get
// ~~~
template<unsigned V, unsigned D, typename T>
inline vector<V, T>& vector_field<V, D, T>::get(const vector<D, unsigned>& pos)
{
  unsigned p=0;
  for(unsigned i=0; i<D; i++)
    p += pos.x[i]*offset[i];

  return f[p];
}


// read
// ~~~~
template<unsigned V, unsigned D, typename T>
void vector_field<V, D, T>::read(std::ifstream& in)
{
  in.read((char*)f, elements*sizeof(T)*V);
}


// write
// ~~~~~
template<unsigned V, unsigned D, typename T>
void vector_field<V, D, T>::write(std::ofstream& out) const
{
  out.write((const char*)f, elements*sizeof(T)*V);
}

// swap
// ~~~~
template<unsigned V, unsigned D, typename T>
void vector_field<V, D, T>::swap(vector_field<V, D, T>& v1,
				 vector_field<V, D, T>& v2)
{
  vector<V, T>* f;

  f = v1.f;
  v1.f = v2.f;
  v2.f = f;
}


#endif 
