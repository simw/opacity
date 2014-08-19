#ifndef QPERPARRAYNEW_H
#define QPERPARRAYNEW_H

#include "../Arrays.h"

using namespace SwArrays;

/**
	@class QperpCalculator
	A class to calculate the necessary sums and products of q and k vectors
	necessary to calculate an opacity expansion. This automates the expansion of the
	dot products inside the C's and B's in GLV for arbitrary order in opacity, rather
	than trying to calculate by hand. All sums of q's from i up to
	the length n are needed, dotted into a sum from j to n. Further, the result
	but with one or more of the q's set to zero (ie a delta function in the integral)
	is needed.
	@author Simon Wicks <simon_wicks@yahoo.com>
*/
class QperpCalculator
{
protected:
  /// n = dimension, order in opacity, number of q's to calculate with
  long n;
  /// There are 2^n different permutations of setting each q to zero, pre calculate this
  long twotothen;

  /// Are we dealing with a correlated medium? 0 = no, 1 = yes
  bool correlated;

  /// Store2D the z value, ie the zeroes
  long zin;

  /// The magnitudes of the q vectors
  MyArray qs;
  /// The polar angles of the q vectors (this is q_perp, ie effectively only 2D vector)
  MyArray ths;
  // 1st index = 0 number, i = row, j = column
  /// Array of vector q_i dot q_j
  array2D qiqj;
  /// Array of q_i.q_j summed from i to n and j to n
  array3D sumqiqj_itonjton;

  /// k vector magnitude - the emitted gluon
  double k;
  /// k vector polar angle 
  double thk;
  /// Array of khat dot q_i (ie khat is unit vector, independent of k magnitude)
  array1D khatqi;
  /// Sum of khat.q_i summed from i to n
  array2D sumkhatqi_iton;
  /// Sum of khat.q_i summed from i to n and j to n
  array3D sumkhatqikhatqj_itonjton;
  /// Sum of k.q_i summed from i to n and j to n - ie depdendent on k
  array3D sumkqikqj_itonjton;

public:

  /// Constructor - takes the wanted dimension, number of collisions
  /// TODO?: templatize this, make n effectively known, thus making loops more efficient?
  /// and whether we're dealing with a correlated medium
  QperpCalculator( unsigned long n_, bool correlated_ );
  ~QperpCalculator();

  /// Set the input array of q-vectors, ths is relative to the fixed direction of k
  void SetQperps( MyArray& qs_, MyArray& ths_ );
  /// Once the q-perps are set, calculate all the further matrices involving q_i.q_j
  void CalcQmatrices();
  void CalcQmatrices2();
  /// Once the q-perps are set, calculate all the further matrices involving khat.q_i
  void CalcKhatMatrices();

  /// Set the magnitude of the k vector (direction is always in fixed direction)
  void SetK( double k_ );

  /// Set the z value - ie which 'zeroes' value we want
  void SetZ( unsigned long z_ );
  /// Retrieve a specific q-vector magnitude
  inline double GetQi( unsigned long i ) const;
  /// Retrieve a specific q-vector polar angle
  inline double GetThi( unsigned long i ) const;
  /// Retrieve the result - sum of q_i.q_j from i to (n+1), j to (n+1) including k.q_i
  double GetSumQiQj( unsigned long i, unsigned long j, unsigned long zeroes ) const;
  /// Same as above, but use the value of zeroes already given by the SetZ function
  double GetSumQiQj( unsigned long i_, unsigned long j_ ) const;

};

double QperpCalculator::GetQi( unsigned long i ) const
{
  return qs[i];
}

double QperpCalculator::GetThi( unsigned long i ) const
{
  return ths[i];
}

#endif
