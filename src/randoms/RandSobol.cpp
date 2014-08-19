#include <cmath>

#include "RandSobol.h"
#include "sobol.H"

namespace SwRandoms
{

RandSobol::RandSobol( unsigned long Dimensionality, unsigned long Seed )
  : RandomBase2(Dimensionality), InitialSeed( Seed )
{

}

RandomBase2* RandSobol::clone() const
{
  return new RandSobol( *this );
}

void RandSobol::GetUniforms( SwArrays::MyArray& variates )
{
  unsigned long dim = GetDimensionality();
  double* sobols = new double[dim];
  
  // Produces a vector of quasirandom numbers
  // Updates RunningSeed
  Sobol::i8_sobol( dim, &RunningSeed, sobols );

  // Move our array of floats into the vector
  for (unsigned long j=0; j<dim; j++)
    variates[j] = sobols[j];

  delete [] sobols;
}

void RandSobol::Skip( unsigned long numberOfPaths)
{
  SwArrays::MyArray tmp( GetDimensionality() );
  for (unsigned long j=0; j<numberOfPaths; j++)
    GetUniforms(tmp);
}

void RandSobol::SetSeed( unsigned long Seed )
{
  //InitialSeed = Seed;
  // Sobol sequence should start from 0
  InitialSeed = 0;
  RunningSeed = InitialSeed;
  srand48( Seed );
}

void RandSobol::Reset()
{
  srand48( InitialSeed );
}

} // End of namespace SwRandoms
