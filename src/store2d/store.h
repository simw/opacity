#ifndef STORE_H
#define STORE_H

#include <vector>

#include "../Arrays.h"
#include "boost/multi_array.hpp"
#include "statisticsmc.h"
#include "../Wrapper.h"

using namespace StatGathering;
using namespace SwArrays;

class Parameters;

/**
  A two dimensional store of results, with an arbitrary length array of StatisticsMC
  derived objects for each point.
  Mechanisms for storage to file and resumption from file are provided.
  File format: section begins with '# Begin Store2D data' and ends with 
  '# End Store2D data'. In between this are all the data point, the number corresponding
  to the Parameters section in the same file (which has already been read).
  Blank lines and lines beginning with # are ignored.

  Dimension1 = k; Dimension2 = x;

  TODO: generalise this class for an arbitrary dimension result.
*/
class Store2D
{
  /// friend to << so that we can overload it to output a Store2D object
  friend std::ostream& operator<< ( std::ostream& out, const Store2D& store );
  friend std::istream& operator>> ( std::istream& in, Store2D& store );

private:
  /// The number of points in dimension 1, excluding the starting point
  long SizeDim1;
  /// The Starting point of dimension 1
  double MinDim1;
  /// The finishing point of dimension 1
  double MaxDim1;  //
  /** The step size (derived from Size,Min,Max) in dimension 1,  
      MaxDim1 = MinDim1 + StepDim1 * SizeDim1
  */
  double StepDim1;

  /// The number of points in dimension 2, excluding the starting point
  long SizeDim2;
  /// The Starting point of dimension 2
  double MinDim2;
  /// The finishing point of dimension 2
  double MaxDim2;  //
  /** The step size (derived from Size,Min,Max) in dimension 2,  
      MaxDim2 = MinDim2 + StepDim2 * SizeDim2
  */
  double StepDim2;

  /// The size of the array of information for each point
  long SizePerPoint;
  /** The array of pointers (wrappers) to Statistics gatherers
  The first dimension is the 2D store combined into 1 dimension
  The second dimension corresponds to a size of SizePerPoint
  */
  boost::multi_array<Wrapper<StatisticsMC>, 2 > stats;

  /// Helper function to turn (IndexDim1,IndexDim2) into a 1D index
  inline long GetIndex( long IndexDim1, long IndexDim2 ) const;

public:
  Store2D();
  /// Constructor, supplying the dimensions of the store grid
  Store2D( long SizeDim1_, long SizeDim2_, long SizePerPoint_ );

  void SetParameters( Parameters& inParams );
  
  /** Get all the store data from a file, and resume from that point
    @return Integer error code: either 0, or negative number
    indicating error
  */
  int ReadFromFile( std::string FileName_ );
  /// Write everything to file
  void WriteToFile( std::string FileName_, bool append );
  /// Resize the store
  int SetSize( long Size1_, long Size2_, long SizePerPoint_ );

  /// Set the limits in dimension 1
  void SetLimitsDim1( double MinDim1_, double MaxDim1_ );
  /// Set the limits in dimension 2
  void SetLimitsDim2( double MinDim2_, double MaxDim2_ );

  /**
    Add a new Monte Carlo point to our results - at (IndexDim1_,IndexDim2_) with an array
    of values given by Values_
  */
  void AddPoint( long IndexDim1_, long IndexDim2_, MyArray& Values_ );

  /// Get the coordinate in dimension 1 corresponding to the given index
  inline double GetCoordDim1( long IndexDim1_ ) const;
  /// Get the coordinate in dimension 2 corresponding to the given index
  inline double GetCoordDim2( long IndexDim2_ ) const;
  /// Get the 2D coordinates corresponding to the given indices
  inline void GetCoords( long IndexDim1_, long IndexDim2_, 
			double& ValDim1, double& ValDim2 ) const;

  inline long GetLengthDim1() const;

  inline long GetLengthDim2() const;

  long GetNumIterations( long opac, long kth, long xth );
};

// 0 <= kth <= sizeK; 0 <= xth <= sizeX
long Store2D::GetIndex( long IndexDim1_, long IndexDim2_ ) const
{
  return ( IndexDim1_*(SizeDim2+1) + IndexDim2_ );
}

double Store2D::GetCoordDim1( long IndexDim1_ ) const
{
  return ( MinDim1 + static_cast<double>(IndexDim1_) * StepDim1 );
}

double Store2D::GetCoordDim2( long IndexDim2_ ) const
{
  return ( MinDim2 + static_cast<double>(IndexDim2_) * StepDim2 );
}

void Store2D::GetCoords( long IndexDim1_, long IndexDim2_,
				 double& ValDim1, double& ValDim2 ) const
{
  ValDim1 = GetCoordDim1( IndexDim1_ );
  ValDim2 = GetCoordDim2( IndexDim2_ );
}

long Store2D::GetLengthDim1() const
{
  return SizeDim1;
}

long Store2D::GetLengthDim2() const
{
  return SizeDim2;
}


#endif
