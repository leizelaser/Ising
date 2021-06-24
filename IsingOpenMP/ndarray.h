/*
  2D container with pointers on one dimension (for neighbor array)
  p[np] links to f[e2*nf] or null (np >= nf)
*/
#ifndef NDARRAY_H
#define NDARRAY_H

#include <iostream>
#include <fstream>
#include <math.h>

template <class T>
class ndarray
{
  public:
  T * f;               // f[elements1*elements2]
  unsigned elements1;  // 
  unsigned elements2;  //
  unsigned pelements1;  //virtualelements1 >= elements1
  T ** p;              // p[virtualelements1]

  ndarray();
  ndarray(unsigned, unsigned);
  ndarray(unsigned, unsigned, unsigned);
  ~ndarray();
  
  // if resize successfully, return true
  // "resize" clears the content of f, while "extend" keeps them 
  bool resize(unsigned, unsigned);
  bool resize(unsigned, unsigned, unsigned);

  void linkptr(unsigned, unsigned);

  // f is relocated. set p manually
  bool resizep(unsigned);
  bool extendf1(unsigned);
  bool extendf2(unsigned);
  // link p to f (the e1 elements from beginning)
  void clearp();
  void autolink();

  // Access to the pointers
  T* operator[](const unsigned);
  void print();
};

// ndarray
// ~~~~~~~
template <class T>
ndarray<T>::ndarray(): elements1(0), elements2(0), pelements1(0), f(NULL), p(NULL)
{
}

// ndarray
// ~~~~~~~
template <class T>
ndarray<T>::ndarray(unsigned e1, unsigned e2): elements1(e1), elements2(e2), pelements1(e1)
{
  f = new T[e1*e2];
  p = new T*[e1];
  for (int i=0; i<e1; i++)
  {
    p[i] = &(f[i*e2]);
  }
}

// ndarray
// ~~~~~~~
template <class T>
ndarray<T>::ndarray(unsigned pe1, unsigned e1, unsigned e2): elements1(e1), elements2(e2), pelements1(pe1)
{
  if (pe1<e1)
  {
    std::cout<<"p's size is smaller than f's (in initialization), force quit."<<std::endl;
    exit(2);
  }
  f = new T[e1*e2];
  p = new T*[pe1];
  for (int i=0; i<pe1; i++)
  {
    p[i] = NULL;
  }
}

// ~ndarray
// ~~~~~~~~
template <class T>
ndarray<T>::~ndarray()
{
  if (f != NULL) delete[] f;
  if (p != NULL) delete[] p;
}

// resize
// ~~~~~~~~
template <class T>
bool ndarray<T>::resize(unsigned e1, unsigned e2)
{
  elements1 = e1; elements2 = e2; pelements1 = e1;
  if (f != NULL) delete[] f;
  if (p != NULL) delete[] p;

  f = new T[e1*e2];
  p = new T*[e1];
  if (f == NULL || p == NULL)
    return false;

  for (int i=0; i<e1; i++)
  {
    p[i] = &(f[i*e2]);
  }
  return true;
}

// resize
// ~~~~~~~~
template <class T>
bool ndarray<T>::resize(unsigned pe1, unsigned e1, unsigned e2)
{
  elements1 = e1; elements2 = e2; pelements1 = pe1;
  if (pe1<e1)
  {
    std::cout<<"p's size is smaller than f's (in resizing), force quit."<<std::endl;
    exit(2);
  }
  if (f != NULL) delete[] f;
  if (p != NULL) delete[] p;

  f = new T[e1*e2];
  p = new T*[pe1];
  if (f == NULL || p == NULL)
    return false;
  for (int i=0; i<pe1; i++)
  {
    p[i] = NULL;
  }
  return true;
}

// linkptr
// ~~~~~~~
template <class T>
inline void ndarray<T>::linkptr(unsigned p1, unsigned f1)
{
  if (p1 >= pelements1|| f1 >= elements1)
  {
    std::cout<<"overflow detected in linking of 2D array, limit = ("
    <<pelements1<<", "<<elements1 
    << ") , input = ("<< p1 << ", " << f1 << ")" << std::endl;
    exit(2);
  }
  p[p1] = &(f[f1*elements2]);
}

template <class T>
bool ndarray<T>::extendf1(unsigned e1)
{
  T * fnew = new T[e1*elements2];
  if (fnew == NULL)
    return false;

  for (int i=0; i<elements1; i++)
    for (int j=0; j<elements2; j++)
      fnew[i*elements2+j] = f[i*elements2+j];
  delete[] f;
  f = fnew;
  elements1 = e1;
  return true;
}

template <class T>
bool ndarray<T>::extendf2(unsigned e2)
{
  T * fnew = new T[elements1*e2];
  if (fnew == NULL)
    return false;

  for (int i=0; i<elements1; i++)
    for (int j=0; j<elements2; j++)
      fnew[i*e2+j] = f[i*elements2+j];

  delete[] f;
  f = fnew;
  elements2 = e2;
  return true;
}

template <class T>
inline bool ndarray<T>::resizep(unsigned pe1)
{
  delete[] p;
  p = new T*[pe1];
  if (p == NULL)
    return false;
  pelements1 = pe1; clearp();
  return true;
}


template <class T>
void ndarray<T>::clearp()
{
  for (int i=0; i<pelements1; i++)
    p[i] = NULL;
}

template <class T>
void ndarray<T>::autolink()
{
  for (int i=0; i<elements1; i++)
    p[i] = &(f[i*elements2]);
}

// []
// ~~
template <class T>
inline T* ndarray<T>::operator[](const unsigned i)
{
  if (i < pelements1) 
    return p[i];
  else
  {
    std::cout<<"overflow detected in 2D array, limit = "<<pelements1 
    << " , i = "<< i << std::endl;
    return NULL;
  }
}

// print
// ~~~~~
template <class T>
inline void ndarray<T>::print()
{
  for (int i=0; i<pelements1; i++)
  {
    std::cout<< i <<": ";
    if (p[i] == NULL)
    {
      std::cout<<"Empty."<<std::endl;
    }
    else
    {
      for (int j=0; j<elements2-1; j++)
      {
        std::cout<<p[i][j]<<", ";
      }
      std::cout<<p[i][elements2-1]<<std::endl;
    }
  }
}

/* A test:
ndarray<int> a(3,4);

for (int i=0; i<3; i++)
    for (int j=0; j<4; j++)
        a[i][j] = i*4+j;
        
a.print();

ndarray<int> a;
a.resize(4,3,4);
a.autolink();
for (int i=0; i<3; i++)
    for (int j=0; j<4; j++)
        a[i][j] = i*4+j;
a.clearp();
a.linkptr(0,0);
a.linkptr(1,1);
a.linkptr(3,2);
a.print();

*/


#endif 
