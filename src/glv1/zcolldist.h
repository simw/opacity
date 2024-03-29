#ifndef ZCOLLDIST_H
#define ZCOLLDIST_H

#include <vector>
#include "../Arrays.h"
#include "Function.h"
#include "../Wrapper.h"

using namespace SwArrays;

class Parameters;

/**
  Generates the z positions, the positions of the collisions along the path length
  These z positions can be correlated with each other
*/
class ZposGenerator
{
private:
  const unsigned long mSize;
  //MyArray mPositions;
  std::vector<double> mPositions;

  double _loverlambda;
  double _weight;

  double _mu;
  double _temp;

  // The reference funtion - the one in the integral
  Wrapper<Function> refFn;
  // The function from which to sample the points
  Wrapper<Function> sampleFn;

public:
  ZposGenerator( Parameters& params_, long opacity_ );

  void FindRandomPositions( MyArray& Randoms );

  void GetTempsMu2s( MyArray& temps, MyArray& mu2s );
  double GetDeltaZi( int i ) const;
  double GetZweight() const;
};

#endif
