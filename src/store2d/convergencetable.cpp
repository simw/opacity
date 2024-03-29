#include <iostream>
#include <cmath>

#include "convergencetable.h"
#include "../constants.h"

namespace StatGathering
{

using namespace SwUtils;

ConvergenceTable::ConvergenceTable( const Wrapper<StatisticsMC>& Inner_ )
 : Inner( Inner_ ), PathsDone( 0 )
{
  // Use a convenient set of defaults
  Increment = 2;
  StoppingPoint = Increment;
  MaxSets = 20;
}

StatisticsMC* ConvergenceTable::clone() const
{
  return new ConvergenceTable( *this );
}

void ConvergenceTable::GetResultsFromInner( std::vector<std::vector<double> >& results_ ) const
{
  // Get the 2D vector of results from the inner
  // For a StatisticsMean, this is just a 1D vector, but allow for arbitrary
  // 2D size
  std::vector<std::vector<double> > thisResult( Inner->GetResultsSoFar() );

  // Iterate over all the rows of the returned result, store them in order
  // at the end of the supplied table, along with the PathsDone
  for (unsigned long i=0; i<thisResult.size(); i++)
  {
    thisResult[i].push_back(PathsDone);
    results_.push_back(thisResult[i]);
  }
}

void ConvergenceTable::AddOneResult( double result )
{
  // Get the work done by the 'Inner' statistics gatherer
  Inner->AddOneResult( result );
  ++PathsDone;

  // If we have hit a power of Increment, we ask the 'Inner' for the results so far,
  // and add it to our array (along with the number of points that contributed to that result)
  if (PathsDone == StoppingPoint)
  {
    StoppingPoint *= Increment;
    GetResultsFromInner( ResultsSoFar );
  }
}

// Add in one set of results to our convergence table
void ConvergenceTable::AddOneSetOfResults( long Number_, 
		std::vector<std::vector<double> > &ResultsSoFar_ )
{
  // If this sends us past the next StoppingPoint, then store current results
  // as long as we haven't just added the results, and we're not at zero
  // (note: this is ResultsSoFar, not ResultsSoFar_)
  if ( PathsDone + Number_ >= StoppingPoint && 
       PathsDone*Increment != StoppingPoint && 
       PathsDone != 0 )
    GetResultsFromInner( ResultsSoFar );

  // Now, pass on this set to the inner statistics gatherer
  Inner->AddOneSetOfResults( Number_, ResultsSoFar_ );
  // Set our total number of points
  PathsDone += Number_;
  StoppingPoint = _FindNextPowerOfTwo( PathsDone );

  GetResultsFromInner( ResultsSoFar );
}

void ConvergenceTable::Reset()
{
  // Pass on the order to the inner
  Inner->Reset();
  // Set our results to zero
  PathsDone = 0;
  StoppingPoint = Increment;
  ResultsSoFar.clear();
}

void ConvergenceTable::SetResultsSoFar( long Number_, 
	std::vector<std::vector<double> > &ResultsSoFar_ )
{
  // Pass on the results to the inner
  Inner->SetResultsSoFar( Number_, ResultsSoFar_ );
  // Set out total number of points
  PathsDone = Number_;
  // Set the next stopping point at which to record to the convergence table
  StoppingPoint = _FindNextPowerOfTwo( PathsDone );
  // Add this current point to our convergence table
  GetResultsFromInner( ResultsSoFar );
}

std::vector<std::vector<double> > ConvergenceTable::GetResultsSoFar() const
{
  // Start with the results given so far
  std::vector<std::vector<double> > tmp( ResultsSoFar );

  // and now add in the current result from the inner object,
  // If we have just hit a power of increment, the result now is identical
  // with that in the last entry in ResultsSoFar, so we don't have
  // to add this extra point.
  if ( (PathsDone*Increment) != StoppingPoint )
  {
    GetResultsFromInner( tmp );
  }
  
  // So now we have a vector of (mean, sd, number) for powers of increment, and the final total
  return tmp;
}

// Allow us to send ConvergenceTable objects to an ostream
// For ease of output to file
std::ostream& operator<<( std::ostream& out, const ConvergenceTable& conTab )
{
  // First, send the Inner to the stream
  out << *( static_cast<const StatisticsMean*> ( conTab.Inner.GetConstPointer() ) );

  // Now, send the convergence table to the stream
  // One entry of this will duplicate the current state of the inner.
  // First, a separator
  out << " ,";

  // Get all the results from the convergence table
  std::vector<std::vector<double> > thisResult( conTab.GetResultsSoFar() );
  // What is the length of the array of results?
  long len = thisResult.size();
  // How many results are we going to write to the file?
  long imax; ( len > conTab.MaxSets ) ? imax = conTab.MaxSets : imax = len;
  
  // Now, print them to the stream
  out << " " << imax;

  long jmax;
  for ( long i=len-imax; i<len; ++i )
  {
    jmax = thisResult[i].size();
    out << " " << jmax;
    for ( long j=0; j<jmax; ++j )
      out << " " << thisResult[i][j];
  }

  return ( out );
}

// Read in elements of the convergence table
std::istream& operator>>( std::istream& in, ConvergenceTable& conTab )
{
  // First, read in the Inner
  in >> *( static_cast<StatisticsMean*> ( conTab.Inner.GetPointer() ) );

  // Now, restore the state of the convergence table
  // First, the separator
  std::string strTmp;
  in >> strTmp;

  std::vector<std::vector<double> > theResults;
  long len;

  // First, get the number of results that the convergence table is holding
  in >> len;
  theResults.resize( len );

  long jmax;
  // Now, iterate through each set
  for ( long i=0; i<len; ++i )
  {
    in >> jmax;
    theResults[i].resize( jmax );
    for ( long j=0; j<jmax; ++j )
      in >> theResults[i][j];
  }

  conTab.ResultsSoFar = theResults;

  return ( in );
}

} // End of StatGathering namespace
