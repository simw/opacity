#ifndef RAND_SOBOL_H
#define RAND_SOBOL_H

#include "Random3.h"

namespace SwRandoms
{

class RandSobol : public RandomBase2
{
public:
  RandSobol( unsigned long Dimensionality, unsigned long Seed = 1);
  
  virtual RandomBase2* clone() const;
  virtual void GetUniforms( SwArrays::MyArray& variates );
  virtual void Skip( unsigned long numberOfPaths );
  virtual void SetSeed( unsigned long Seed );
  virtual void Reset();
  
private:
  long long int InitialSeed;
  long long int RunningSeed;
};

} // End of namespace SwRandoms

#endif
