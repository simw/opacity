#ifndef RAND_DRAND48_H
#define RAND_DRAND48_H

#include "Random3.h"

namespace SwRandoms
{

class RandDrand48 : public RandomBase2
{
public:
  RandDrand48( unsigned long Dimensionality, unsigned long Seed = 1);
  
  virtual RandomBase2* clone() const;
  virtual void GetUniforms( SwArrays::MyArray& variates );
  virtual void Skip( unsigned long numberOfPaths );
  virtual void SetSeed( unsigned long Seed );
  virtual void Reset();
  
private:
  unsigned long InitialSeed;
};

} // End of namespace SwRandoms

#endif
