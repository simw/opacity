//
// C++ Interface: qperpgenerator
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef QPERPGENERATOR3_H
#define QPERPGENERATOR3_H

#include <boost/array.hpp>
#include "../sw_templatemeta.h"
#include "../constants.h"

// Unnamed namespace
namespace
{

const double pi = 3.141592653589793238;

}

/**
  @author 
*/
template<std::size_t n>
class QperpGenerator3
{
private:
  double qmax;
  double mu2;
  boost::array<double, TPower<2,n>::value> weights;
public:
  QperpGenerator3();
  void SetParameters( double temp, double energy, double mu );
  //void GetQsThetas( boost::array<double, 3*n>& randoms, boost::array<double, n>& Qs,
  void GetQsThetas( boost::array<double, 2*n>& randoms, boost::array<double, n>& Qs,
    boost::array<double, n>& Thetas );
  double GetQsEventWeight( long zeroes ) const;
};

template<std::size_t n>
QperpGenerator3<n>::QperpGenerator3()
{

}

template<std::size_t n>
void QperpGenerator3<n>::SetParameters( double temp, double energy, double mu )
{
  qmax = sqrt( 6. * temp * energy );
  mu2 = mu*mu;
}

template<std::size_t n>
//void QperpGenerator3<n>::GetQsThetas( boost::array<double, 3*n>& randoms, 
void QperpGenerator3<n>::GetQsThetas( boost::array<double, 2*n>& randoms, 
  boost::array<double, n>& Qs, boost::array<double, n>& Thetas )
{
  for ( std::size_t i=0; i<n; ++i )
  {
    //Qs[i] = randoms[n+i]*qmax;
    //Thetas[i] = randoms[2*n+i]*2.*pi;
    Qs[i] = randoms[i]*qmax;
    Thetas[i] = randoms[n+i]*2.*pi;
  }

  boost::array<double, n> qweights;
  for ( std::size_t i=0; i<n; ++i )
  {
    qweights[i] = 2.*Qs[i] * mu2 / pow( Qs[i]*Qs[i] + mu2, 2 );
    qweights[i] *= qmax;
  }

  std::vector<bool> _isZeroed(n);
  for ( long j=0; j< TPower<2,n>::value; ++j )
  {
    weights[j] = 1.;
    SwUtils::_NumberToBoolArray( j, _isZeroed, n );
    for ( std::size_t i=0; i<n; ++i )
    {
      if ( !_isZeroed[i] )
        weights[j] *= qweights[i];
      else
        weights[j] *= -1.;
    }
  }
}

template<std::size_t n>
double QperpGenerator3<n>::GetQsEventWeight( long zeroes ) const
{
  return weights[zeroes];
}

#endif
