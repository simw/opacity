//
// C++ Interface: qperpcalculator
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef QPERPCALCULATOR3_H
#define QPERPCALCULATOR3_H

#include "../sw_templatemeta.h"
#include "../constants.h"
#include <boost/array.hpp>

template<std::size_t n>
void _Convert2DPolarToCartesian( boost::array<double, n> &qs, 
  boost::array<double, n> &thetas, boost::array<double, n> &qxs, boost::array<double, n> &qys )
{
  for ( std::size_t i=0; i<n; ++i )
  {
    qxs[i] = qs[i] * cos( thetas[i] );
    qys[i] = qs[i] * sin( thetas[i] );
  }
}

/**
  @author 
*/
template<std::size_t n>
class QperpCalculator3
{
private:
  /// Array of Sum_(i=1 to n, j = k to n)_Qi Dot Qj
  boost::array<boost::array<double, n>, TPower<2,n>::value> SumQs_1i;
  boost::array<boost::array<double, n>, TPower<2,n>::value> SumQs_ii;
  boost::array<boost::array<double, n>, TPower<2,n>::value> SumKhatQ_i;

  boost::array<boost::array<double, n+1>, TPower<2,n>::value> SumQsIncK_1i;
  boost::array<boost::array<double, n+1>, TPower<2,n>::value> SumQsIncK_ii;

public:
  QperpCalculator3();

  void SetK( double k );
  void SetQsThetas( boost::array<double, n> Qs, boost::array<double, n> Thetas );

  double GetSumQskk( long m, long zeroes ) const
  {
    return ( SumQsIncK_ii[zeroes][m-1] );
  };

  double GetSumQs1k( long m, long zeroes ) const
  {
    return ( SumQsIncK_1i[zeroes][m-1] );
  };

};

template<std::size_t n>
QperpCalculator3<n>::QperpCalculator3()
{

}

template<std::size_t n>
void QperpCalculator3<n>::SetK( double k )
{
  const double k2 = k*k;
  boost::array<double, n> sumkq_i;
  for ( long z=0; z<TPower<2,n>::value; ++z )
  {
    for ( std::size_t i=0; i<n; ++i )
    {
      sumkq_i[i] = k * SumKhatQ_i[z][i];
      SumQsIncK_1i[z][i] = SumQs_1i[z][i] + sumkq_i[0] + sumkq_i[i] + k2;
      SumQsIncK_ii[z][i] = SumQs_ii[z][i] + 2.*sumkq_i[i] + k2;
    }
    SumQsIncK_1i[z][n] = sumkq_i[0] + k2;
    SumQsIncK_ii[z][n] = k2;
  }
}

template<std::size_t n>
void QperpCalculator3<n>::SetQsThetas( boost::array<double, n> Qs, 
  boost::array<double, n> Thetas )
{
  boost::array<double, n> qxs, qys;
  _Convert2DPolarToCartesian<n>( Qs, Thetas, qxs, qys );

  boost::array<double, n> sumqxs, sumqys;
  bool notZeroed;

  for ( std::size_t z=0; z<TPower<2,n>::value; ++z )
  {
    notZeroed = !SwUtils::_TestBitI( z, 0 );
    sumqxs[n-1] = notZeroed*qxs[n-1]; sumqys[n-1] = notZeroed*qys[n-1];
    for ( std::size_t i=1; i<n; ++i )
    {
      std::size_t ind = n-1-i;
      notZeroed = !SwUtils::_TestBitI( z, i );
      sumqxs[ind] = sumqxs[ind+1] + notZeroed*qxs[ind];
      sumqys[ind] = sumqys[ind+1] + notZeroed*qys[ind];
    }

    for ( std::size_t i=0; i<n; ++i )
    {
      SumQs_1i[z][i] = sumqxs[0]*sumqxs[i] + sumqys[0]*sumqys[i];
      SumQs_ii[z][i] = sumqxs[i]*sumqxs[i] + sumqys[i]*sumqys[i];
      SumKhatQ_i[z][i] = -sumqxs[i];
    }
  }
}

#endif
