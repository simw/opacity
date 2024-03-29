#ifndef RANDOM3_H
#define RANDOM3_H

#include "../Arrays.h"

namespace SwRandoms
{

class RandomBase2
{
public:
  RandomBase2( unsigned long Dimensionality_ ) : Dimensionality( Dimensionality_ ) {};
  
  inline unsigned long GetDimensionality() const;
  
  virtual RandomBase2* clone() const = 0;
  virtual void GetUniforms( SwArrays::MyArray& variates ) = 0;
  virtual void Skip( unsigned long numberOfPaths ) = 0;
  virtual void SetSeed( unsigned long Seed ) = 0;
  virtual void Reset() = 0;
  
private:
  unsigned long Dimensionality;
};

unsigned long RandomBase2::GetDimensionality() const
{
  return Dimensionality;
}

} // End of SwRandoms namespace

#endif
