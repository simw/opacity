
#ifndef QPERPSTORE2_H
#define QPERPSTORE2_H

#include "boost/multi_array.hpp"

/**
	@author 
*/
class QperpStore2
{
public:
  QperpStore2( long sizeXY_, long sizeZ_ );

  inline void SetZ( long zin_ );
  void SetVal( long x_, long y_, double val_ );
  double GetVal( long x_, long y_ );

private:
  long setZ;
  boost::multi_array<double, 3> vals;
};

QperpStore2::QperpStore2( long sizeXY_, long sizeZ_ )
  : vals( boost::extents[sizeXY_][sizeXY_][sizeZ_] ), setZ( 0 )
{

}

void SetZ( long zin_ )
{
  setZ = zin_;
}

void SetVal( long x_, long y_, long val_ )
{
  vals[x_][y_][setZ] = val_;
}

double GetVal( long x_, long y_ )
{
  return vals[x_][y_][setZ];
}


#endif
