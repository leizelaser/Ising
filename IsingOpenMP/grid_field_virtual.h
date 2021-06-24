#ifndef GRID_FIELD_VIRTUAL_H
#define GRID_FIELD_VIRTUAL_H

#include "vector.h"


// ======================================================================
// grid_field_virtual
// ======================================================================

// This class makes it easier to do VecToInd and IndToVec
// No real data contained
// It can be initialized with a grid_field
template<unsigned D>
class grid_field_virtual {

 public:
  vector<D, unsigned> size;           // number of grid points for each dimension
  vector<D, unsigned> offset;
 public:
  unsigned elements;
 
 public:
  grid_field_virtual();
  grid_field_virtual(const vector<D, unsigned>&);
  grid_field_virtual(const unsigned);

  // make a virtual copy of grid_field: set_size(g.size); 
  // grid_field_virtual(grid_field<D, T>&);
  ~grid_field_virtual();

  void set_size(const vector<D, unsigned>&);
  void set_size(const unsigned);

  //inverse indexing:
  vector<D, unsigned> IndtoVec(unsigned);
  unsigned VectoInd (const vector<D, unsigned>&);
};


// grid_field_virtual
// ~~~~~~~~~~~~~~~~~~
template<unsigned D>
grid_field_virtual<D>::grid_field_virtual():elements(0)
{
}


// grid_field_virtual
// ~~~~~~~~~~~~~~~~~~
template<unsigned D>
grid_field_virtual<D>::grid_field_virtual(const vector<D, unsigned>& s)
{
  set_size(s);
}

template<unsigned D>
grid_field_virtual<D>::grid_field_virtual(const unsigned s)
{
  set_size(s);
}

// ~grid_field_virtual
// ~~~~~~~~~~~~~~~~~~~
template<unsigned D>
grid_field_virtual<D>::~grid_field_virtual()
{;}


// set_size
// ~~~~~~~~
template<unsigned D>
void grid_field_virtual<D>::set_size(const vector<D, unsigned>& s)
{
  size = s;
  elements = 1;
  for(unsigned i=0; i<D; i++) {
    offset[i] = elements;
    elements *= size.x[i];
  }
}


// set_size
// ~~~~~~~~
template<unsigned D>
void grid_field_virtual<D>::set_size(const unsigned s)
{
  vector<D, unsigned> square;

  for(unsigned k=0; k<D; k++)
    square[k] = s;
  set_size(square);
}


// indexing
// ~~~~~~~~

template<unsigned D>
inline vector<D, unsigned> grid_field_virtual<D>::IndtoVec(unsigned ind)
{
  vector<D, unsigned> vec;
  for (unsigned i=D-1;i>0;i--)
  {
    vec[i] = ind/offset[i];
    ind = ind - vec[i] * offset[i];
  }
  vec[0] = ind/offset[0];
  return vec;
}

template<unsigned D>
inline unsigned grid_field_virtual<D>::VectoInd(const vector<D, unsigned>& pos)
{
  unsigned p=0;
  for(unsigned i=0; i<D; i++)
    p += pos.x[i]*offset[i];
  return p;
}

#endif
