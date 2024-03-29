//
// C++ Interface: qperpcalculator1
//
// Description: A wrapper to use the old version of QperpCalculator in the new GlvRadiative3
//
// Author: Simon Wicks <simonw@phys.columbia.edu>
//
#ifndef QPERPCALCULATOR1_H
#define QPERPCALCULATOR1_H

#include <boost/array.hpp>

#include "../sw_templatemeta.h"
#include "../Arrays.h"
#include "../glv1/qperparraynew.h"

/**
  @author Simon Wicks <simonw@phys.columbia.edu>
  A wrapper around the QperpCalculator from the old version, version1, to use
  in the newer GlvRadiative3 template calcaluation
  The interface necessary for GlvRadiative3 is converted into the relevant calls
  for the old code. In this way, the old, tested qperpcalculator can be used in the
  newer code, to compare to results for any newer versions.
*/
template<std::size_t n>
class QperpCalculator1
{
private:
  QperpCalculator qperps;
public:
  QperpCalculator1() : qperps( n, true )
  {

  };

  void SetK( double k )
  {
    qperps.SetK( k );
  };

  void SetQsThetas( boost::array<double, n>& Qs, 
    boost::array<double, n>& Thetas )
  {
    SwArrays::MyArray qs1(n), thetas1(n);
    for ( std::size_t i=0; i<n; ++i )
    {
      qs1[i] = Qs[i]; thetas1[i] = Thetas[i];
    }
    qperps.SetQperps( qs1, thetas1 );
  };

  double GetSumQskk( long m, long zeroes ) const
  {
    return qperps.GetSumQiQj( m, m, zeroes );
  };

  double GetSumQs1k( long m, long zeroes ) const
  {
    return qperps.GetSumQiQj( 1, m, zeroes );
  };

};

#endif
