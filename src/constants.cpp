#include <vector>
#include "constants.h"

namespace SwUtils
{

long _FindNextPowerOfTwo( long value )
{
  long myreturn = 0;
  --value;
  do
  {
    ++value;
    if ( _IsPowerOfTwo( value ) )
      myreturn = value;
  }
  while( myreturn == 0 );

  return myreturn;
}

long _Combinatoric( long n, long r )
{
  // This returns the combinatoric factor n C r
  // ie n! / ( (n-r)! r! )
  // Note: sum of n C r over r gives 2^n

  if ( r >= n || r <= 0 )
    return 1;

  // _factorial (n,b) = n!/(b-1)! ie n*(n-1)*...*b
  return ( _factorial( n, n-r+1 ) / _factorial( r ) );
}

} // End of SwUtils namespace
