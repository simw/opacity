#include <cmath>
#include <iostream>
#include <algorithm>

#include "zcolldist.h"
#include "../parameters.h"
#include "../constants.h"

using namespace SwUtils;

ZposGenerator::ZposGenerator( Parameters& params_, long opacity_ )
  : mSize( opacity_ ), mPositions( opacity_+1 )
{
  std::vector<double> ReturnedParamsDouble;
  std::vector<long> ReturnedParamsLong;
  std::list<std::string> ReturnedParamsString;

  // First, we get the parameters of the medium
  // This is: 1st = mu, 2nd = temperature, 3rd = gluon mass, 4th = gluon mean free path
  // Here, we only need the gluon mass
  // the zDist will need the others
  ReturnedParamsDouble = params_.GetParametersDouble( "@mediumParams" );
  _mu = ReturnedParamsDouble[0];
  _temp = ReturnedParamsDouble[1];
  double _gluonlambda = ReturnedParamsDouble[3];

  // The jet path length in the medium
  ReturnedParamsDouble = params_.GetParametersDouble( "@pathLength" );
  double _length = ReturnedParamsDouble[0];

  double _maxlen;
  _maxlen = _length * 5.;
  _loverlambda = _length / _gluonlambda;

  MyArray zmin(1), zmax(1), param(1);

  // Reference function => the function under the integral
  // Sample function => the function from which to sample our Monte Carlo points

  // For uniform reference function
//  UniformFunction refFn1( 1 );
//  param[0] = 0.; zmin[0] = 0.; zmax[0] = _length;
//  refFn1.SetLimits( zmin, zmax); refFn1.SetParams( param );

  // For exponential decay reference function
  ExpDecayFunction refFn1( 1 );
  param[0] = _length / static_cast<double>(mSize+1); zmin[0] = 0.; zmax[0] = _maxlen;
  refFn1.SetLimits( zmin, zmax ); refFn1.SetParams( param );

  // For uniform sample function
//  UniformFunction sampFn1( 1 );
//  param[0] = 0.; 
//  sampFn1.SetParams( param );

  // For exponential decay sample function
  ExpDecayFunction sampFn1( 1 );
  param[0] = _length / static_cast<double>(mSize+1);
  sampFn1.SetParams( param );


  refFn = refFn1;
  sampleFn = sampFn1;

  mPositions[0] = 0.;
}

void ZposGenerator::FindRandomPositions( MyArray& Randoms )
{
  MyArray zmin(1), zmax(1), zIn(1), zOut(1), refOut(1), sampOut(1);
  //_weight = _loverlambda;
  _weight = 1. / _factorial( mSize );
  for (unsigned long i=0; i<mSize; ++i)
  {
    //zmin[0] = mPositions[i]; zmax[0] = refFn->GetXmax(0);
    zmin[0] = 0.; zmax[0] = refFn->GetXmax(0);
    sampleFn->SetLimits( zmin, zmax );
    zIn[0] = Randoms[i]; sampleFn->NormedIntegralInverse( zIn, zOut ); 
    mPositions[i+1] = mPositions[i] + zOut[0]; //mPositions[i+1] = zOut[0];
    refFn->function( zOut, refOut ); sampleFn->function( zOut, sampOut );
    _weight *= _loverlambda;
    _weight *= refOut[0] / sampOut[0];
  }

}

void ZposGenerator::GetTempsMu2s( MyArray& temps, MyArray& mu2s )
{
  for (unsigned long i=0; i<mSize; ++i)
  {
    temps[i] = _temp; mu2s[i] = _mu*_mu;
  }
}

double ZposGenerator::GetDeltaZi( int i ) const
{
  double tmp;
  tmp = mPositions[i] - mPositions[i-1];
  
  return tmp;
}

double ZposGenerator::GetZweight() const
{
  return _weight;
}

