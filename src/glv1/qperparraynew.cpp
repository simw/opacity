
#include "qperparraynew.h"
#include "../constants.h"

using namespace SwUtils;

/// A number of utility functions, in an unnamed namespace - ie for use only in this file
namespace
{

// FillRow functions: these start from the matrix qiqjs and calculates a matrix
// that is the sums of qiqjs - the (i,j)th element is the sum of all qiqjs
// in the box (i->n, j->n). In the interference calculation, these are the only
// sums that are necessary to calculate.
// There are many different ways of doing this, but here is an implementation 
// that only evaluates 0.5*n(n+1) sums.

// Plus: we need to calculate the same matrix, but with specific qi's set to zero
// For a correlated medium (eg variable density), there are 2^n different permutations
// For an uncorrelated geometry (fixed density), we don't need to calculate anything new
// for each permutation - just keep the original, but set the last qi's to zero

// _FillFirstRow: fills the (n-1)th row and column of sumqiqjs
void _FillFirstRow( unsigned long n, unsigned long z, std::vector<bool> zeroed,
				 		array2D& qiqjs, array3D& sumqiqjs )
{
  // If the (n-1)th row is zeroed, jump ahead to the answer
  if (zeroed[n-1])
  {
    for (unsigned long i=0; i<n; ++i)
    {
      sumqiqjs[z][n-1][i] = 0.; sumqiqjs[z][i][n-1] = 0.;
    }
  }
  else
  {
    // The bottom right hand corner of the matrix has just one contribution
    sumqiqjs[z][n-1][n-1] = qiqjs[n-1][n-1];

    // Now fill in the lower row and the right most column
    for (unsigned long i=2; i<(n+1); ++i)
    {
      sumqiqjs[z][n-1][n-i] = sumqiqjs[z][n-1][n-i+1] + (!zeroed[n-i])*qiqjs[n-1][n-i];
      sumqiqjs[z][n-i][n-1] = sumqiqjs[z][n-1][n-i];
    }
  }
}

// _FillNextRow - having done the summations for the bottom row and the right column,
// spread out the results to the next row and column
void _FillNextRow( unsigned long n, unsigned long i, 
	unsigned long z, std::vector<bool> zeroed, array2D& qiqjs, array3D& sumqiqjs )
{
  for (unsigned long j=0; j<(i+1); ++j)
    sumqiqjs[z][i][i-j] = sumqiqjs[z][i+1][i-j] 
	+ sumqiqjs[z][i][i-j+1] - sumqiqjs[z][i+1][i-j+1] 
	+ (!zeroed[i])*(!zeroed[i-j])*qiqjs[i][i-j];

  for (unsigned long j=1; j<(i+1); ++j)
    sumqiqjs[z][i-j][i] = sumqiqjs[z][i][i-j];
}

// _FillAllRows - driving function for FillFirstRow and FillNextRow
// We could probably do this by a recursive function that calls itself, but
// let's do it the simple way.
void _FillAllRows2DSymmetric( unsigned long n, array2D& qiqjs, array3D& sumqiqjs, long z )
{
  std::vector<bool> zeroed( n );
  _NumberToBoolArray( z, zeroed, n ); 
  _FillFirstRow( n, z, zeroed, qiqjs, sumqiqjs );
  for (unsigned long i=2; i<(n+1); ++i)
    _FillNextRow( n, n-i, z, zeroed, qiqjs, sumqiqjs );
}

void _FillAllRows1D( unsigned long n, array1D& khatqis, array2D& sumkhatqis, long z )
{
  // We will be stepping through setting each qi to zero
  // Create a vector of bools to indicate which qi is zero for this iteration
  std::vector<bool> zeroed( n );

  // Convert number z to an array of bools, the ith element saying whether qi should be zero
  _NumberToBoolArray( z, zeroed, n ); 
  // Starting point for the iteration: (n-1)th has just one element, khat.q_{n-1}
  sumkhatqis[z][n-1] = (!zeroed[n-1])*khatqis[n-1];
  // Now iterate over the rest, adding one khatqi at each step
  for (unsigned long i=2; i<(n+1); ++i)
    sumkhatqis[z][n-i] = sumkhatqis[z][n-i+1] + (!zeroed[n-i])*khatqis[n-i];
}

void _Fill2DFrom1D( unsigned long n, array2D& sumkhatqi_iton, 
			array3D& sumkhatqikhatqj_itonjton, long z )
{
  for (unsigned long i=0; i<n; ++i)
  {
    for (unsigned long j=0; j<n; ++j)
    {
      sumkhatqikhatqj_itonjton[z][i][j] = sumkhatqi_iton[z][i] + sumkhatqi_iton[z][j];
    }
    sumkhatqikhatqj_itonjton[z][i][n] = sumkhatqi_iton[z][i];
  }
  for (unsigned long j=0; j<n; ++j)
  {
    sumkhatqikhatqj_itonjton[z][n][j] = sumkhatqi_iton[z][j];
  }
  sumkhatqikhatqj_itonjton[z][n][n] = 0.;
}

// _DotProducts2D: want qi.qj. Take n vectors, expressed in polar coords
// ( magnitudes mags, polar vectors ths ), and fill the 2 dimension nxn matrix
// res with the results qi.qj
void _DotProducts2D( unsigned long n, MyArray& mags, MyArray& ths, array2D& res )
{
  // TODO: put in asserts to check that the dimensions of mags, ths and res are correct

  // The operations sin and cos are numerically expensive
  // For the dot products: will need cos( ths[i] - ths[j] )
  // = cos(ths[i]) cos(ths[j]) + sin(ths[i]) sin(ths[j])
  // Here: evaluate them 2n times, instead of possibly n^2 or 0.5*n(n+1)
  MyArray sines(n), cosines(n);
  for (unsigned long i=0; i!=n; ++i)
  {
    sines[i] = sin( ths[i] );
    cosines[i] = cos( ths[i] );
  }

  // Now, find the 2D matrix qi.qj
  // We know that it is symmetric, so first we'll evaluate the off-diagonal on
  // one side, then mirror it over to the other side
  for (unsigned long i=0; i!=n; ++i)
    for (unsigned long j=i+1; j!=n; ++j)
    {
      res[i][j] = mags[i] * mags[j] * ( cosines[i]*cosines[j] + sines[i]*sines[j] );
      res[j][i] = res[i][j];
    }

  // Now all that's left is the diagonal
  for (unsigned long i=0; i!=n; ++i)
    res[i][i] = mags[i] * mags[i];

  // And we're done!
}


// _DotProducts1D: want k.qi. Take n vectors, expressed in polar coords
// ( magnitudes mags, polar vectors ths ) and dot into one fixed vector (mag1, th1),
// to give the 1 dimensional array res
void _DotProducts1D( unsigned long n, double mag1, double th1, 
			MyArray& mags, MyArray& ths, array1D& res )
{
  // TODO: put in asserts to check that the dimensions of mags, ths and res are correct

  // No need to pre-calculate the cosine this time, we'll just jump in
  for (unsigned long i=0; i!=n; ++i)
    res[i] = mag1 * mags[i] * cos( th1 - ths[i] );

  // And we're done!
}


}
// End of the unnamed namespace

QperpCalculator::QperpCalculator( unsigned long n_, bool correlated_ )
  : n(n_), correlated(correlated_), 
    qs( n_ ), ths( n_ ),
    qiqj( boost::extents[n_][n_] ),
    thk(pi),
    khatqi( boost::extents[n_] )
{
  twotothen = power(2,n);

  if ( correlated==0 )
  {
    sumqiqj_itonjton.resize( boost::extents[1][n+1][n+1] );
    sumkhatqi_iton.resize( boost::extents[1][n] );
    sumkhatqikhatqj_itonjton.resize( boost::extents[1][n+1][n+1] );
    sumkqikqj_itonjton.resize( boost::extents[1][n+1][n+1] );
  }
  else
  {
    sumqiqj_itonjton.resize( boost::extents[twotothen][n+1][n+1] );
    sumkhatqi_iton.resize( boost::extents[twotothen][n] );
    sumkhatqikhatqj_itonjton.resize( boost::extents[twotothen][n+1][n+1] );
    sumkqikqj_itonjton.resize( boost::extents[twotothen][n+1][n+1] );
  }
}

QperpCalculator::~QperpCalculator()
{

}

void QperpCalculator::SetZ( unsigned long z_ )
{
  zin = z_;

  if ( correlated )
  {

  }

}

double QperpCalculator::GetSumQiQj( unsigned long i_, unsigned long j_, unsigned long z_ ) const
{
  long i = i_-1; long j = j_-1; long z = z_;

  if ( !correlated )
  {
    if (i < z) i = z_;
    if (j < z) j = z_;
    z = 0;
  }

  boost::array<array3D::index,3> idx = {{z,i,j}};
  return ( sumqiqj_itonjton( idx ) + sumkqikqj_itonjton( idx ) );
}

void QperpCalculator::SetK( double k_ )
{
  k = k_; double k2 = k*k;

  long zmax;
  correlated ? zmax = twotothen-1 : zmax = 0;
  
  for( index3D z=0; z<=zmax; ++z )
  {
    for ( index3D i=0; i<n+1; ++i )
    {
      for ( index3D j=0; j<n+1; ++j )
      {
        sumkqikqj_itonjton[z][i][j] = sumkhatqikhatqj_itonjton[z][i][j] * k + k2;
      }
    }
  }
}

void QperpCalculator::SetQperps( MyArray& qs_, MyArray& ths_ )
{
  qs = qs_;
  ths = ths_;

  CalcQmatrices2();
}

// CalcQmatrices2: once we have received a set of q vectors, calculate all the qiqjs
// and sumqiqjs, as well as khatqis and sumkhatqis
void QperpCalculator::CalcQmatrices2()
{
  // Here, we pre-calculate all the necessary qi.qj elements and their useful sums
  // We'll also calculate the angular part of k.qi, but we'll keep the magnitude of k
  // generic at the moment. This way, we can iterate over all possible k magnitudes
  // for one set of qi.qj's.

  // First, calculate all qi.qj elements. This fills 2D matrix qiqj with qi.qj 
  _DotProducts2D( n, qs, ths, qiqj );
  // Now calculate k.qi elements for unit vector k. This fills khatqi.
  _DotProducts1D( n, 1., thk, qs, ths, khatqi );

  // Will iterate over all sets of zeroed out q's, if we have a correlated medium
  long zmax;
  correlated ? zmax = twotothen-1 : zmax = 0;

  for (long z=0; z<=zmax; ++z)
  {
    // Now we want to fill the 2D matrix sumqiqj_itonjton.
    _FillAllRows2DSymmetric( n, qiqj, sumqiqj_itonjton, z );
    // Finally, filling sumkhatqi_iton matrix
    _FillAllRows1D( n, khatqi, sumkhatqi_iton, z );
    // Convert the 1D khatqi_iton into a 2D khatqikhatqj ready for adding to qiqj_itonjton
    _Fill2DFrom1D( n, sumkhatqi_iton, sumkhatqikhatqj_itonjton, z );
  }
}

// The old way

namespace
{

void _OuterProduct1D1D( unsigned long n, MyArray& mags1, MyArray& ths1, MyArray& mags2, MyArray& ths2, array2D& res )
{
  for (unsigned long i=0; i!=n; ++i)
    for (unsigned long j=i+1; j!=n; ++j)
    {
      res[i][j] = mags1[i] * mags2[j] * cos( ths1[i] - ths2[j] );
      res[j][i] = res[i][j];
    }
  
  for (unsigned long i=0; i!=n; ++i)
    res[i][i] = mags1[i] * mags2[i];
}

void _OuterProduct0D1D( unsigned long n, double mag1, double th1, MyArray& mags2, MyArray& ths2, array1D& res )
{
  for (unsigned long i=0; i!=n; ++i)
    res[i] = mag1 * mags2[i] * cos( th1 - ths2[i] );
}



void _Iterate2DForZeroes( unsigned long n, array2D& qiqjs, array3D& sumqiqjs )
{
  // Now we go through all the possible combinations of zeroed columns and rows
  std::vector<bool> zeroed( n );
  double tmp;
  for (unsigned long z=1; z<pow(2,n); ++z)
  {
    // First, calculate the array of zeroed bools
    _NumberToBoolArray( z, zeroed, n );
    // Now, iterated over all the elements
    for (unsigned long i=0; i<n; ++i) 
    { 
      for (unsigned long j=0; j<n; ++j) 
      {
        // If either the row or column has been zeroed, subtract off from answer
        tmp = sumqiqjs[0][i][j];
        for (unsigned long ii=i; ii<n; ii++)
        {
          for (unsigned long jj=j; jj<n; jj++)
          {
            if ( zeroed[ii] || zeroed[jj] )
              tmp -= qiqjs[ii][jj];
          }
        }
        sumqiqjs[z][i][j] = tmp;
      }
    }
  }
}

void _Iterate1DForZeroes( unsigned long n, array1D& khatqis, array2D& sumkhatqis )
{
  std::vector<bool> zeroed(n);
  double tmp;
  for (unsigned long z=1; z<pow(2,n); ++z)
  {
    // First, calculate the array of zeroed bools
    _NumberToBoolArray( z, zeroed, n );
    // Now, iterated over all the elements
    for (unsigned long i=0; i<n; ++i) 
    { 
      // Now deal with k part
      tmp = sumkhatqis[0][i];
      for (unsigned long ii=i; ii<n; ++ii)
      {
        if ( zeroed[ii] )
        {
          tmp -= khatqis[ii];
        }
      }
      sumkhatqis[z][i] = tmp;
    }
  }
}

}

void QperpCalculator::CalcQmatrices()
{
  // Array of q_i.q_j summed from i to n
  array2D sumqiqj_iton;
  // Array of q_i.q_j summed from j to n
  array2D sumqiqj_jton;

  sumqiqj_iton.resize( boost::extents[n][n] );
  sumqiqj_jton.resize( boost::extents[n][n] );

  // Here, we pre-calculate all the necessary qi.qj elements and their useful sums

  // First, calculate all qi.qj elements and khat.qi
  _OuterProduct1D1D( n, qs, ths, qs, ths, qiqj );
  _OuterProduct0D1D( n, 1., thk, qs, ths, khatqi );
  
  // Now, calculate sumqiqj_jton[i][j] = Sum[ qiqj[i][jj], {jj, j, n} ];
  // ie sumqiqj_jton[i][j] = sumqiqj_iton[i][j+1] + qiqj[i][j]
  // sumqiqj_iton[i][j] = transpose sumqiqj_jton[i][j] = sumqiqj_jton[j][i]
  for (index2D i=0; i!=n; ++i) 
  {
    sumqiqj_jton[i][n-1] = qiqj[i][n-1];
    sumqiqj_iton[n-1][i] = sumqiqj_jton[i][n-1];
    for (index2D j=(n-2); j>=0; --j)
    {
      sumqiqj_jton[i][j] = sumqiqj_jton[i][j+1] + qiqj[i][j];
      sumqiqj_iton[j][i] = sumqiqj_jton[i][j];
    }
  }
  
  // Now we calculate the final product
  
  // sumqiqj_itonjton[i][j] = Sum[ sumqiq_jton[ii][j], {ii, i, n} ];
  // ie sumqiqj_itonjton[i][j] = sumqiqj_itonjton[i+1][j] + sumqiqj_jton[i][j]
  sumqiqj_itonjton[0][n][n] = 0.;
  for (long j=0; j<n; ++j) 
  {
    sumqiqj_itonjton[0][n][j] = 0.;
    sumqiqj_itonjton[0][n-1][j] = sumqiqj_jton[n-1][j];
    for (long i=(n-2); i>=0; --i)
    {
      sumqiqj_itonjton[0][i][j] = sumqiqj_itonjton[0][i+1][j] + sumqiqj_jton[i][j];
    }
  }
  
  // Now throw in sumkhatqi_iton[i] = Sum[ khatqi[ii], {ii,i,n} ];
  // ie sumkhatqi_iton[i] = sumkhatqi_iton[i+1] + khatqi[i]
  sumkhatqi_iton[0][n-1] = khatqi[n-1];
  for (long i=(n-2); i>=0; --i)
  {
    sumkhatqi_iton[0][i] = sumkhatqi_iton[0][i+1] + khatqi[i];
  }

  if ( correlated )
  {
    // Now we go through all the possible combinations of zeroed columns and rows
    _Iterate2DForZeroes( n, qiqj, sumqiqj_itonjton );
    _Iterate1DForZeroes( n, khatqi, sumkhatqi_iton );
  }

  // Now calc sumkhatqikhatqj_itonjton
  long zmax; correlated ? zmax = twotothen-1 : zmax = 0;
  for (long z=0; z<=zmax; ++z)
  {
    for (long i=0; i<n; ++i)
    {
      for (long j=0; j<n; ++j)
      {
        sumkhatqikhatqj_itonjton[z][i][j] = sumkhatqi_iton[z][i] + sumkhatqi_iton[z][j];
      }
      sumkhatqikhatqj_itonjton[z][i][n] = sumkhatqi_iton[z][i];
    }
    for (long j=0; j<n; ++j)
    {
      sumkhatqikhatqj_itonjton[z][n][j] = sumkhatqi_iton[z][j];
    }
    sumkhatqikhatqj_itonjton[z][n][n] = 0.;
  }
}

