#ifndef QPERPDIST_H
#define QPERPDIST_H

#include <vector>
#include <valarray>
#include "../Arrays.h"
#include "qperparraynew.h"
#include "Function.h"
#include "../Wrapper.h"

using namespace SwArrays;

// This class generates randomly distributed qs and thetas
// from any desired distribution, sends the generated qs and thetas
// into the QperpCalculator object, and keeps track of the
// weight relative to a fixed Gyulassy-Wang distribution

class QperpGenerator
{
private:
  /// Number of dimensions of q, theta
  const unsigned long _dim;
  /// The class dealing with the details of the calculation, once filled with qs, thetas
  QperpCalculator qps;
  /// The weight of the event compared to the fixed Gyulassy-Wang distribution
  MyArray _weights;
  /// The product of all the weights, ie the event weight
  double _totWeight;

  /// Whether the medium is correlated
  bool _correlated;

  /// The code for which qs are zeroed - int converts to binary gives which qs are zeroed
  unsigned long _zeroSet;
  /// The total number of zeroed qs
  unsigned long _numZeroedQs;
  /// The array of zeroed qs - ie zeroSet in binary form, with 0<->1 inverted
  std::vector<bool> _isZeroed;
  /// Whether we have an even number of zeroed qs (hence a -1 in the event weight)
  bool _evenZeroes;

  /// The distribution function from which to choose the qs
  Wrapper<Function> qDistFn;
  /// The reference function in the integral - yukawa (ie Gyulassy-Wang)
  Wrapper<Function> yukDistFn;
  /// The distribution function from which to choose the thetas
  Wrapper<Function> thDistFn;
  /// The reference function in the integral - uniform
  Wrapper<Function> uniDistFn;

public:
  /// Constructor: feed in the number of dimensions, and whether the medium is correlated
  QperpGenerator( unsigned long n_, bool correlated_ );
  /// Destructor
  ~QperpGenerator();

  /// Supply random numbers in inForQs, inForThs, generates random qs between qmins and qmaxs, ths between thmins and thmaxs - using the Debye masses from mu2s
  void FindRandomQs( MyArray& inForQs, MyArray& qmins, MyArray& qmaxs, 
                           MyArray& inForThs, MyArray& thmins, MyArray& thmaxs, 
                           MyArray& mu2s );
  /// Set the combination of zeroed qs from an integer
  void SetZeroedQs( unsigned long num );
  /// Set the k value
  inline void SetK( double k );

  /// Get out the SumQiQj - passed straight onto QperpCalculator object
  inline double GetSumQiQj( unsigned long i, unsigned long j ) const;
  /// Get the weight of the event, ie the ratio of the distribtutions from which the qs, ths are generated and the ones in the integrals
  double GetQeventWeight( ) const;
  /// Ask whether a specific q is zeroed
  inline bool IsZeroed( unsigned long i );

};

bool QperpGenerator::IsZeroed( unsigned long i )
{
  return _isZeroed[i];
}

void QperpGenerator::SetK( double k )
{
  qps.SetK( k );
}

double QperpGenerator::GetSumQiQj( unsigned long i, unsigned long j ) const
{
  return qps.GetSumQiQj( i, j, _zeroSet );
}

#endif
