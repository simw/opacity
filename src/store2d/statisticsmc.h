#ifndef STATISTICSMC_H
#define STATISTICSMC_H

#include <vector>

namespace StatGathering
{

/**
  Base class for statistics gatherer for a Monte Carlo integrator
  based on code from Joshi 'C++ Design Patterns and Derivatives Pricing'

  This base class is pure virtual, defines an interface to interact with the statistics gatherer
  This represents one 'point': each point is filled with a single double, but can return
  a 2D vector with information about its results.
*/
class StatisticsMC
{
public:
  StatisticsMC() {}
  virtual ~StatisticsMC() {}

  /// clone method, to interface with the 'Wrapper' template
  virtual StatisticsMC* clone() const = 0;

  /** Add one result to the mix
      This is the standard way of adding in another result
  */
  virtual void AddOneResult( double result ) = 0;
  /** Add many results to the mix
    This might be from merging two sets of results or similar
    \param Number The total number of results to add in
    \param ResultsSoFar_ A 2D vector of doubles, which gives the relevant structure of the results
         (this structure is not specified in this base class)
  */
  virtual void AddOneSetOfResults( long Number, 
	std::vector<std::vector<double> > &ResultsSoFar_ ) = 0;

  /// Reset all internal values, start again
  virtual void Reset() = 0;
  /** Set the results as supplied
    \param Number The total number of results so far
    \param ResultsSoFar_ A 2D vector of doubles, which gives the relevant structure of the results
         (this structure is not specified in this base class)
  */
  virtual void SetResultsSoFar( long Number_, 
	std::vector<std::vector<double> > &ResultsSoFar_ ) = 0;
  /** Get a suitably structured 2D vector of the results so far
      This base class does not specify the structure to return
  */
  virtual std::vector<std::vector<double> > GetResultsSoFar() const = 0;

private:

};


/**
  A simple implementation of a statistics gatherer
  Keeps a running tab on the total sum of results and the sum of results^2
  The results: 1D vector - the mean, the standard deviation
*/
class StatisticsMean : public StatisticsMC
{
  /// friend to << so that we can overload the << operator to output a ConvergenceTable object
  friend std::ostream& operator<< ( std::ostream& out, const StatisticsMean& stats );
  friend std::istream& operator>> ( std::istream& in, StatisticsMean& stats );

public:
  StatisticsMean();

  virtual StatisticsMC* clone() const;

  virtual void AddOneResult( double result );
  /**
    Add in a large set of results to what we already have
    \param Number_ The number of results to add in
    \param ResultsSoFar_ 2D vector of results size (1,2)
	, with two elements, [0] = mean, [1] = standard deviation
  */
  virtual void AddOneSetOfResults( long Number_, 
	std::vector<std::vector<double> > &ResultsSoFar_ );

  /// Reset all internal values to zero
  virtual void Reset();
  /// Set to a specific set of results
  virtual void SetResultsSoFar( long Number_, 
	std::vector<std::vector<double> > &ResultsSoFar_ );
  /**
    Get the results so far
    2D vector of results: in this case, dimensions (1,2)
    The two elements are: mean, standard deviation
    \return Return by value, a 2D vector of dimension (1,2) giving mean, standard deviation
  */
  virtual std::vector<std::vector<double> > GetResultsSoFar() const;

private:
  /// The sum of all the input results
  double RunningSum;
  /// The sum of all the squares of input results
  double RunningSum2;
  /// The total number of points input
  unsigned long PathsDone;
};

} // End of StatGathering namespace

#endif
