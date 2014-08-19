//
// C++ Interface: qperpstore1
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef QPERPSTORE1_H
#define QPERPSTORE1_H

/**
	@author 
*/
class QperpStore1
{
public:
  QperpStore1( long sizeXY_, long sizeZ_ );
  ~QperpStore1();

  inline void SetZ( long zin_ );
  void SetVal1( long x_, long y_, double val_ );
  void SetVal2( long x_, long y_, double val_ );
  void MakeVals3( 
  double GetVal( long x_, long y_ );

private:
  double* vals1;
  double* vals2;
  double* vals3;

  long sizeXY;
  long sizeXY2;
  long sizeZ;

  long bookmarkZ;

  inline long XYtoIndex( long x_, long y_ );
  void Reset();
};

QperpStore1::QperpStore1( long sizeXY_, long sizeZ_ )
{
  sizeXY = sizeXY_; 
  sizeXY2 = sizeXY*sizeXY;
  sizeZ = sizeZ_;

  vals = 0;
  long totSize = sizeZ * sizeXY2;
  vals1 = new double( totSize );
  vals2 = new double( totSize );
  vals3 = new double( totSize );

  bookmarkZ = 0;
}

QperpStore1::~QperpStore1()
{
  delete [] vals1;
  vals1 = 0;
  delete [] vals2;
  vals2 = 0;
  delete [] vals3;
  vals3 = 0;
}

void SetZ( long zin_ )
{
  bookmarkZ = sizeXY2 * zin_;
}

void SetVal( long x_, long y_, long val_ )
{
  vals[ XYtoIndex( x_, y_ ) ] = val_;
}

double GetVal( long x_, long y_ )
{
  return vals[ XYtoIndex( x_, y_ ) ];
}

long XYtoIndex( long x_, long y_ )
{
  return ( bookmarkZ + y_*sizeXY + x_ );
}

void Reset()
{
  long imax = sizeZ * sizeXY2;
  
  for (long i=0; i<imax; ++i)
    vals[i] = 0.;
}

#endif
