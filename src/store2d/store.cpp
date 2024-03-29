#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>

#include "store.h"
#include "convergencetable.h"
#include "statisticsmc.h"
#include "../parameters.h"
#include "../constants.h"

Store2D::Store2D()
{

}

Store2D::Store2D( long SizeDim1_, long SizeDim2_, long SizePerPoint_ )
{
  SetSize( SizeDim1_, SizeDim2_, SizePerPoint_ );
}

void Store2D::SetParameters( Parameters& MyParameters )
{
  std::vector<double> ReturnedParamsDouble1, ReturnedParamsDouble2;
  std::list<std::string> ReturnedParamsString;  

  // First, get the details from the parameters object
  double min1, max1, min2, max2;
  long size1, size2;
  ReturnedParamsDouble1 = MyParameters.GetParametersDouble( "@storedim1" );
  if ( ReturnedParamsDouble1.size() < 3 )
    std::cerr << "Error with @storedim1" << std::endl;
  ReturnedParamsDouble2 = MyParameters.GetParametersDouble( "@storedim2" );
  if ( ReturnedParamsDouble2.size() < 3 )
    std::cerr << "Error with @storedim2" << std::endl;
  size1 = long(ReturnedParamsDouble1.at(2)+0.5);
  size2 = long(ReturnedParamsDouble2[2]+0.5);
  min1 = ReturnedParamsDouble1[0];
  min2 = ReturnedParamsDouble2[0];
  max1 = ReturnedParamsDouble1[1];
  max2 = ReturnedParamsDouble2[1];

  SetSize( size1, size2, 1 );
  SetLimitsDim1( ReturnedParamsDouble1[0], ReturnedParamsDouble1[1] );
  SetLimitsDim2( ReturnedParamsDouble2[0], ReturnedParamsDouble2[1] );
}

int Store2D::SetSize( long Size1_, long Size2_, long SizePerPoint_ )
{
  SizeDim1 = Size1_;
  SizeDim2 = Size2_;
  SizePerPoint = SizePerPoint_;

  stats.resize( boost::extents[ (SizeDim1+1)*(SizeDim2+1) ][SizePerPoint] );

  // We want to set up all the StatisticsMC as pointers to
  // ConvergenceTable objects containing StatisticsMean objects
  StatisticsMean gatherer;
  ConvergenceTable tab(gatherer);

  for (long i=0; i!=(SizeDim1+1)*(SizeDim2+1); ++i)
    for (long opac=0; opac!=SizePerPoint; ++opac)
      stats[i][opac] = tab;

  return 0;
}

void Store2D::SetLimitsDim1( double MinDim1_, double MaxDim1_ )
{
  MinDim1 = MinDim1_;
  MaxDim1 = MaxDim1_;
  // For convenience, we calculate the fixed step size and store it
  StepDim1 = (MaxDim1-MinDim1)/double(SizeDim1);
}

void Store2D::SetLimitsDim2( double MinDim2_, double MaxDim2_ )
{
  MinDim2 = MinDim2_;
  MaxDim2 = MaxDim2_;
  // For convenience, we calculate the fixed step size and store it
  StepDim2 = (MaxDim2_-MinDim2_)/static_cast<double>(SizeDim2);
}

void Store2D::AddPoint( long IndexDim1_, long IndexDim2_, MyArray& Values_ )
{
  for (long i=0; i<SizePerPoint; ++i)
    stats[GetIndex( IndexDim1_, IndexDim2_ )][i]->AddOneResult( Values_[i] );
}

long Store2D::GetNumIterations( long opac, long kth, long xth )
{
  std::vector<std::vector<double> > res = 
	stats[GetIndex( kth, xth )][opac]->GetResultsSoFar();
  
  std::vector<std::vector<double> >::iterator it;
  it = res.end(); --it;
  return static_cast<long>( it->at(2) );
}

int Store2D::ReadFromFile( std::string FileName_ )
{
  // Reading in a file, and putting it in the store
  // The code here is very similar to that in the parameters.cpp for reading
  // in the parameters from a file

  // First, we open the file and check that it exists
  std::ifstream FileIn;
  FileIn.open( FileName_.c_str(), std::ios::in );
  if ( FileIn.fail() )
  {
    std::cerr << "Store2D, unable to open file " << FileName_;
    std::cerr << " for reading." << std::endl;
    return -20;
  }

  // Each line, we'll read into LineReadIn
  // Each line of settings will be counted in NumberOfLines
  std::string LineReadIn;

  // Ok, we have an open file, we want to find the data section
  // The string that identifies the beginning of the settings section
  std::string BeginString = "# Begin Store2D data";
  // The string that identifies the end of the settings section
  std::string EndString = "# End Store2D data";

  // Indicator whether we have found the settings section yet
  bool FoundStart = false;
  // Run through the file until we find the settings section
  do
  {
    if ( FileIn.eof() )
    {
      std::cerr << "Store2D, reached end of file at line before "; 
      std::cerr << "finding data section" << std::endl;
      return -21;
    }
    getline( FileIn, LineReadIn );
    boost::trim( LineReadIn );
    if ( LineReadIn == BeginString )
      FoundStart = true;
  }
  while ( !FoundStart );

  // Now, pass off the logic to the >> operator
  FileIn >> *this;

  // Now we're done, we close the file
  FileIn.close();
  
  // Everything has worked fine. 
  // Return a success code
  return 0;
}

void Store2D::WriteToFile( std::string FileName_, bool append )
{
  // First, open the file
  std::ofstream FileOut;
  if ( append )
  {
    FileOut.open( FileName_.c_str(), std::ios::app );
  }
  else
  {
    FileOut.open( FileName_.c_str(), std::ios::trunc );
  }
  
  // Check for errors
  if ( FileOut.fail() )
  {
    std::cerr << "Store2D, unable to open file " << FileName_;
    std::cerr << " for writing." << std::endl;
  }

  FileOut << "# Begin Store2D data" << std::endl;
  FileOut << *this;
  FileOut << "# End Store2D data" << std::endl;
  FileOut.close();
}

std::ostream& operator<<( std::ostream& out, const Store2D& store )
{
  // Assume that all the preliminaries have been written
  // All we have to do here is write the data to file / stream

  // Format: dim1 dim2 results....
  // Incrementing dim1, then dim2
  
  // If there are multiple data sets per point, output one whole set then the next
  // Hence, outer iteration is over data sets per point
  for (long i=0; i<store.SizePerPoint; ++i)
  {
    // Inner iteration over all points
    for (long dim2=0; dim2<store.SizeDim2; ++dim2)
    {
      for (long dim1=0; dim1<store.SizeDim1; ++dim1)
      {
        // First, write the k, x coords to file
        // (for ease of reading either manually or by eg Mathematica)
        out << store.GetCoordDim1( dim1 ) << " " << store.GetCoordDim2( dim2 ) << " ";
        // Pass off the logic of writing to the ConvergenceTable object
        const ConvergenceTable* conTab = static_cast<const ConvergenceTable*>
                 ( store.stats[ store.GetIndex( dim1, dim2 ) ][i].GetConstPointer() );
        out << *conTab;
        out << std::endl;
      }
    }
  }

  // We're done, all overloading of << has to return the supplied ostream
  return ( out );
}

std::istream& operator>>( std::istream& in, Store2D& store )
{
  // TODO: improve this to handle possible alterations to the file
  // eg blank lines, comment lines etc, handle errors and bad files properly

  double coordDim1, coordDim2;
  // Hence, outer iteration is over data sets per point
  for (long i=0; i<store.SizePerPoint; ++i)
  {
    // Inner iteration over all points
    for (long dim2=0; dim2<store.SizeDim2; ++dim2)
    {
      for (long dim1=0; dim1<store.SizeDim1; ++dim1)
      {
        // First, we have the coords in dim1 and dim2
        in >> coordDim1 >> coordDim2;
        // Pass off the logic of reading to the ConvergenceTable object
        ConvergenceTable* conTab = static_cast<ConvergenceTable*>
                 ( store.stats[ store.GetIndex( dim1, dim2 ) ][i].GetPointer() );
        in >> *conTab;
      }
    }
  }

  return ( in );
}

