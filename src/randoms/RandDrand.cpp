#include <cmath>

#include "RandDrand.h"

namespace SwRandoms
{

RandDrand48::RandDrand48( unsigned long Dimensionality, unsigned long Seed )
  : RandomBase2(Dimensionality), InitialSeed( Seed )
{

}

RandomBase2* RandDrand48::clone() const
{
  return new RandDrand48( *this );
}

void RandDrand48::GetUniforms( SwArrays::MyArray& variates )
{
  for (unsigned long j=0; j<GetDimensionality(); j++)
    variates[j] = drand48();
}

void RandDrand48::Skip( unsigned long numberOfPaths)
{
  SwArrays::MyArray tmp( GetDimensionality() );
  for (unsigned long j=0; j<numberOfPaths; j++)
    GetUniforms(tmp);
}

void RandDrand48::SetSeed( unsigned long Seed )
{
  InitialSeed = Seed;
  srand48( Seed );
}

void RandDrand48::Reset()
{
  srand48( InitialSeed );
}

} // End of namespace SwRandoms
