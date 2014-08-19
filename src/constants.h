#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <vector>

namespace SwUtils
{

const double pi = 3.141592653589793238;
const double overpi = 0.3183098861837906715;
const double hbarc = 0.197;

// _factorial (n,b) = n!/(b-1)! ie n*(n-1)*...*b
inline long _factorial( long value, long endval = 1 )
{
  if (value > endval)
    return ( value * _factorial( value - 1 ) );

  return endval;
}

long _FindNextPowerOfTwo( long value );
long _Combinatoric( long n, long r );

inline long power( long num, long pow )
{
  long res = num;
  for (long i=1; i!=pow; ++i)
    res *= num;

  return res;
}

inline bool _TestBitI( unsigned long num, unsigned long i )
{
  int bit = ((num >> i) & 1);
  if (bit == 1)
    return true;
    
  return false;
}

inline void _NumberToBoolArray( unsigned long num, std::vector<bool>& boolarray, unsigned long length )
{
  for (unsigned long i=0; i!=length; ++i)
  {
    boolarray[(length-1)-i] = _TestBitI( num, i );
  }
}

inline long _CountNumZeroes( unsigned long num, unsigned long length )
{
  long tot = 0;
  for (unsigned long i=0; i!=length; ++i)
  {
    if ( _TestBitI( num, i ) )
      ++tot;
  }
  return tot;
}

inline bool _IsPowerOfTwo ( long value )
{
  if (value < 1)
    return false;

  return (value & (~value+1)) == value;  //~value+1 equals a two's complement -value
}

} // End of SwUtils namespace

#endif

