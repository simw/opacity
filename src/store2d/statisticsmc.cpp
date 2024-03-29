#include <iostream>
#include <cmath>
#include "statisticsmc.h"

namespace StatGathering
{

StatisticsMean::StatisticsMean()
  : RunningSum(0.), RunningSum2(0.), PathsDone(0)
{

}

void StatisticsMean::AddOneResult( double result )
{
  // Adding in one result
  // Hence, increment PathsDone by one, add the result to the running sum,
  // add the squre to RunningSum2
  PathsDone++;
  RunningSum += result;
  RunningSum2 += result*result;
}

void StatisticsMean::AddOneSetOfResults( long Number_, 
	std::vector<std::vector<double> > &ResultsSoFar_ )
{
  PathsDone += Number_;

  // ResultsSoFar_ is passed in as two elements: the mean and the standard deviation
  double mean, sd, num;
  num = static_cast<double>(Number_);
  mean = ResultsSoFar_[0][0];
  sd = ResultsSoFar_[0][1];

  // ... but we store internally as total total of squares
  // now convert from mean, sd to total, total of squres
  RunningSum += mean * num;
  RunningSum2 += (sd*sd*num + mean*mean) * num;
}

void StatisticsMean::Reset()
{
  PathsDone = 0;
  RunningSum = 0.;
  RunningSum2 = 0.;
}

void StatisticsMean::SetResultsSoFar( long Number_, 
	std::vector<std::vector<double> > &ResultsSoFar_ )
{
  Reset();
  AddOneSetOfResults( Number_, ResultsSoFar_ );
}

std::vector<std::vector<double> > StatisticsMean::GetResultsSoFar() const
{
  // Construct our 2D vector
  std::vector<std::vector<double> > Results(1);
  Results[0].resize(2);

  // Calculate the mean and standard deviation
  Results[0][0] = RunningSum / PathsDone;
  Results[0][1] = sqrt( (RunningSum2 / PathsDone - Results[0][0]*Results[0][0]) / PathsDone );

  return Results;
}

StatisticsMC* StatisticsMean::clone() const
{
  return new StatisticsMean(*this);
}

std::ostream& operator<< ( std::ostream& out, const StatisticsMean& stats )
{
  std::vector<std::vector<double> > theResults = stats.GetResultsSoFar();
  double mean = theResults[0][0];
  double sd = theResults[0][1];

  out << mean << " " << sd << " " << stats.PathsDone;

  return ( out );
}

std::istream& operator>> ( std::istream& in, StatisticsMean& stats )
{
  std::vector<std::vector<double> > theResults( 1 );
  theResults[0].resize(2);
  unsigned long pathsDone;

  double mean;
  double sd;

  in >> mean >> sd >> pathsDone;
  theResults[0][0] = mean;
  theResults[0][1] = sd;

  stats.SetResultsSoFar( pathsDone, theResults );
 
  return ( in );
}

} // End of StatGathering namespace
