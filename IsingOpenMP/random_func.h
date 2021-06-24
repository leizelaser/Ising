/*
  Random library based on rand() and RAND_MAX
*/
#ifndef RANDOM_FUNC_H
#define RANDOM_FUNC_H

#include <iostream>
#include <fstream>
#include <math.h>

// uniform distribution between [0,1]
inline double rand_unif()
{
  double result;
  while (true)
  {
    result = double(rand())/RAND_MAX;
    if (result < 1)
      return result;
  }
}

// uniform integer in 0, 1, ... n-1
inline unsigned rand_int(unsigned n)
{
  //return unsigned(rand_unif()*n);
  return rand()%n;
}


#ifndef VECTOR_H
//#warning "Include 'vector.h' before 'random_func.h' for better performance"
#else
// randvec
// ~~~~~~~
template <unsigned D>
inline vector<D, double> RandDoubleVec()
{
  vector<D, double> result;
  for (int i=0; i<D; i++)
    result[i] = rand_unif();
  return result;
} 
#endif

// choose n random number from 0 to Nlimit-1, fill them to a given address
// when n<<Nlimit
/* bool RandFill(unsigned Nlimit, unsigned n, unsigned *result)
{
  if (Nlimit<n)
    {
      std::cout<<"Error in RandFill: Nlimit("<<Nlimit<<") >= n("<< n <<")."<<std::endl;
      return false;
    }  
  unsigned temp[n]; //sorted (from small to large) result

  result[0] = rand_int(Nlimit);
  result[1] = rand_int(Nlimit-1);
  if (result[1]>=result[0])
    result[1] += 1;

  if (n>2)
  {
    for(int i=2; i<n; i++)
    {
      bool flag = true;
      while(flag)
      {
        result[i] = rand_int(Nlimit);
        for (int j=0; j<i; j++)
          if (result[i] == result[j])
          {
            flag = true;
            continue;
          }
        if (flag == false) break;
      }
    }
  }
  return true;
}
*/

// choose 2 random numbers from 0 to Nlimit-1, fill them to a given address
bool RandFillTwo(unsigned Nlimit, unsigned *result)
{
  if (Nlimit<2)
    {
      std::cout<<"Error in RandFill: Nlimit("<<Nlimit<<") >= 2."<<std::endl;
      return false;
    }  
  result[0] = rand_int(Nlimit); result[1] = rand_int(Nlimit-1);
  if (result[1]>=result[0]) result[1] += 1;
  return true;
}

// choose 2 random numbers from 0 to Nlimit-1, fill them to a given address
bool Shuffle(unsigned Nlimit, unsigned *result)
{
  unsigned swapind, t;
  for (int i=0; i<Nlimit; i++)
    result[i] = i;
  int i=Nlimit-1;
  while (i>1)
  {
    swapind = rand_int(Nlimit);
    if (swapind == i)
      continue;
    else
    {
      t = result[i];
      result[i] = result[swapind]; 
      result[swapind] = t;
    }
    i--;
  }
  return true;
}

#endif 
