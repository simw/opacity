#ifndef CONVERGENCETABLE_H
#define CONVERGENCETABLE_H

#include "statisticsmc.h"
#include "../Wrapper.h"

namespace StatGathering
{

/**
  An implementation of a statistics gatherer, that keeps track of the results given in powers of (Increment)
  It provides the final mean, standard deviation, and also the mean and sd at any power of (Increment)
  below the final result
*/
class ConvergenceTable : public StatisticsMC
{
  /// friend to << so that we can overload the << operator to output a ConvergenceTable object
  friend std::ostream& operator<< ( std::ostream& out, const ConvergenceTable& conTab );
  friend std::istream& operator>> ( std::istream& in, ConvergenceTable& conTab );

public:
  ConvergenceTable( const Wrapper<StatisticsMC>& Inner_ );

  virtual StatisticsMC* clone() const;

  /**
    Add one result. Passes along that result to the inner, if we have hit 'StoppingPoint' then
    ask for the mean and sd from inner, and store in our array
  */
  virtual void AddOneResult( double result );
  /// Add in a set of results
  virtual void AddOneSetOfResults( long Number_, 
	std::vector<std::vector<double> > &ResultsSoFar_ );

  /// Reset to an empty table
  virtual void Reset();
  /// Pass in the results of a simulation so far, one final set
  virtual void SetResultsSoFar( long Number_, 
	std::vector<std::vector<double> > &ResultsSoFar_ );
  /** Get a 2D vector of all our results so far
    \return 2D vector, formatted as for Inner, for all the powers of increment passed, 
    and one for where we are now, in ascending order of number
  */
  virtual std::vector<std::vector<double> > GetResultsSoFar() const;

private:
  /**
    The 'Inner' does the calculating from point to point, then passes up the information to 
    the convergence table when asked - in the function GetResultsFromInner
  */
  Wrapper<StatisticsMC> Inner;
  /** Our 2D vector of results so far for each power of increment that we have passed, 
      where number is the number of points that contributed (ie increment^n).
      The specific format of this vector depends on the type of StatisticsMC used.
  */
  std::vector<std::vector<double> > ResultsSoFar;
  /// The next point at which we need to ask for the results from the inner, and then store them
  long StoppingPoint;
  /// Number of points evaluated so far (also stored in the 'Inner')
  long PathsDone;
  /// The power at which we want to store results (eg at 2,4,8,16 ... or 10,100,1000 ...)
  long Increment;
  /// The maximum number of sets to record to file at the end
  long MaxSets;
  
  /** Ask the inner for a set of results, add it to the end of results_  */
  void GetResultsFromInner( std::vector<std::vector<double> >& results_ ) const;
};

/**
  Outputs the table to the stream in the format: "Mean SD Number"
  for each of the powers of Increment passed so far, in reverse order
  (ie the final answer comes first, works down to smaller numbers)
*/
std::ostream& operator<<( std::ostream& out, const ConvergenceTable& conTab );
std::istream& operator>>( std::istream& in, const ConvergenceTable& conTab );

}

#endif
