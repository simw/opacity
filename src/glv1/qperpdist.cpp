#include <cmath>
#include <iostream>

#include "qperpdist.h"
#include "Function.h"
#include "../constants.h"

using namespace SwUtils;

QperpGenerator::QperpGenerator( unsigned long n, bool correlated_ )
  : _dim( n ), qps( n, correlated_ ), _weights( n ), _isZeroed( n+1 )
{
  YukawaFunction yuk( _dim );
  UniformFunction uni( _dim );
  
  qDistFn = yuk;
  yukDistFn = yuk;
  thDistFn = uni; 
  uniDistFn = uni;

  _correlated = correlated_;

  // k is never zeroed
  _isZeroed[n] = false;
}

QperpGenerator::~QperpGenerator()
{

}

void QperpGenerator::SetZeroedQs( unsigned long num )
{
  _zeroSet = num;

  if ( _correlated )
  {
    _NumberToBoolArray( num, _isZeroed, _dim );

    _numZeroedQs = 0; _evenZeroes = true;
    for (unsigned long i=0; i<_dim; ++i) 
    {
      if (_isZeroed[i])
      {   ++_numZeroedQs; _evenZeroes = !_evenZeroes; }
    }
  }
  else
  {
    for (unsigned long i=0; i<num; ++i)
      _isZeroed[i] = 1;
    for (unsigned long i=num; i<_dim; ++i)
      _isZeroed[i] = 0;

    _numZeroedQs = num;
    _evenZeroes = !(num % 2);
  }
}

void QperpGenerator::FindRandomQs( MyArray& inForQs, MyArray& qmins, MyArray& qmaxs, 
                           MyArray& inForThs, MyArray& thmins, MyArray& thmaxs, 
                           MyArray& mu2s )
{
  // We have inputs: inForQs, inForThs - randomly distributed numbers between 0,1
  // We want to produce Qs, Ths distributed between qmins,qmaxs and thmins,thmaxs
  MyArray _qs( _dim ), _ths( _dim );
  SetZeroedQs( 0 );

  qDistFn->SetLimits( qmins, qmaxs );
  yukDistFn->SetLimits( qmins, qmaxs );
  thDistFn->SetLimits( thmins, thmaxs );
  uniDistFn->SetLimits( thmins, thmaxs );

  qDistFn->SetParams( mu2s ); 
  yukDistFn->SetParams( mu2s );

  qDistFn->NormedIntegralInverse( inForQs, _qs );
  thDistFn->NormedIntegralInverse( inForThs, _ths );

  qps.SetQperps( _qs, _ths );

  MyArray tmp1( _dim ), tmp2( _dim ), tmp3( _dim ), tmp4( _dim );
  qDistFn->function( _qs, tmp1 );
  yukDistFn->function( _qs, tmp2 );
  thDistFn->function( _ths, tmp3 );
  uniDistFn->function( _ths, tmp4 );

  _totWeight = 1.;
  for (unsigned long i=0; i<_dim; ++i)
  {
    _weights[i] = tmp2[i]/tmp1[i] * tmp4[i]/tmp3[i];
    _totWeight *= _weights[i];
  }
}

double QperpGenerator::GetQeventWeight() const
{
  double tmp = 1.;

  for (unsigned long i=0; i<_dim; i++)
  {
    if ( !_isZeroed[i] )
      tmp *= _weights[i];
  }
  if ( !_evenZeroes ) tmp = -tmp;
  return tmp;
}
