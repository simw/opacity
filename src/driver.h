//
// C++ Interface: driver
// Description: A generic driver for a Monte Carlo calculation
//
// Author: Simon Wicks <simon_wicks@yahoo.com>
//
#ifndef DRIVER_H
#define DRIVER_H

#include <string>
#include <vector>

#include "Wrapper.h"
#include "randoms/Random3.h"
#include "randoms/RandDrand.h"
#include "randoms/RandSobol.h"
#include "Arrays.h"
#include "parameters.h"

/**
  @author Simon Wicks <simonw@phys.columbia.edu>
  A generic driver for a Monte Carlo calculation
  At run time, different random number generator can be used: for example
  pseudorandoms such as standard drand48(); quasirandoms such as Sobol sequences (the latter
  can vastly increase the rate of convergence of Monte Carlo integration).

  Template parameters: TcalcEngine, Tstore, numOfRandoms
  TcalcEngine is an object that performs the calculation, with the interface
    .SetParameters( Parameters inParams )
    .SetRandoms( boost::array<numOfRandoms> randomNumbers )
    .SetCoord( long coord, double value )
    .GetAnswer( std::vector<double> answers )
  Tstore is an object in which to keep the results, with the interface
    .SetParameters( Parameters inParams )
    .AddPoint( long dim1, long dim2, std::vector<double> result )
    .WriteToFile( std::string filename, bool appendToFile )
  numOfRandoms is the number of random numbers to be used in the simulation

  TODO: add the number of dimensions of the result as a template parameter,
        generalize the iteration over deterministic dimensions to be n-dimensional
        add this number of deterministic dimensions as a template parameter

  TODO: add multithreading, so that multiple calculation engines can be run in parallel
        while adding into the same data store
*/
template<typename TcalcEngine, typename Tstore, std::size_t numOfRandoms>
class Driver
{
private:
  /// Supplied template parameter type, to perform the calculation
  TcalcEngine calculator;
  /// Supplied template parameter type, to store the results of the calculation
  Tstore storage;
  /// Object to read runtime parameters from file, supplied to calculation and storage
  Parameters myParameters;
  /// Random number generator, on the interval [0,1)
  Wrapper<SwRandoms::RandomBase2> randomGenerator;
  
public:
  Driver();
  ~Driver();
  /** Setup routine - read in runtime parameters, set up all objects as necessary
    @resume ="y" to resume from a previous run, ="n" to start a new one
    @inputFile Filename of a) previous run file if resuming (including runtime parameters
    on the end), or b) runtime parameters file if not resuming
  **/
  int Setup( std::string resume, std::string inputFile );
  /// Read any runtime parameters for the Driver object from the Parameters object
  void SetParameters( Parameters &myParameters );
  /// Run a single iteration of the monte carlo simulation
  void RunOneIteration();
  /// Save the stored results to a file
  void SaveResults( std::string outputFile );
};

template<typename TcalcEngine, typename Tstore, std::size_t numOfRandoms>
Driver<TcalcEngine, Tstore, numOfRandoms>::Driver()
  : myParameters( "Opacity3" )
{
  // Nothing to do on creation, most setup done by the Setup() function
  // so that errors on setup can be understood

  // The only thing is the constructor of myParameters,
  // to tell it that the version of this program is Opacity version 3
  // TODO: make the 'Opacity3' section describe all necessary parts of the program
  // for a compatible read in / out
}

template<typename TcalcEngine, typename Tstore, std::size_t numOfRandoms>
Driver<TcalcEngine, Tstore, numOfRandoms>::~Driver()
{
  // Nothing to do in the destructor
}

template<typename TcalcEngine, typename Tstore, std::size_t numOfRandoms>
int Driver<TcalcEngine, Tstore, numOfRandoms>::Setup( std::string resume, std::string inputFile )
{
  // Function to do all necessary run time setup of objects
  // TODO: give a meaningful return code, so if this fails the program quits

  // The runtime settings are read in from a file with a specific format
  // This format is governed by the 'Parameters' object, see there for specifics

  // First, read in the file into the object
  // If we are resuming from a previous run, the parameters have been written to the
  // end of the file, and the Parameters objects can find them there
  // All the tokens in the file are read into a 'map'
  myParameters.ParseInputFile( inputFile );

  // Pass the parameters object to each object in turn
  // First, this one, the driver
  SetParameters( myParameters );
  // Second, the storage object - eg for size of each dimensions, number of points etc 
  storage.SetParameters( myParameters );
  // Third, the calculation object - eg GLV specific params, mu, temperature etc
  calculator.SetParameters( myParameters );

  // If we are resuming from a previous run, we also need to read in the statistics
  // from that run into the data store
  if ( resume == "y" )
    storage.ReadFromFile( inputFile );

  // Finally, check that all the inputs from the parameters file were used at some point
  // Otherwise, the user might think that some parameters are being set that really aren't
  int allUsed = myParameters.CheckForUnaccessedParameters();
  if ( allUsed != 0 )
  {
    std::cerr << "Driver::Setup - " << allUsed << " unused parameters from input file.\n";
    return 1;
  }
  else
  {
    std::cout << "All parameters read successfully from " << inputFile << std::endl;
  }

  // TODO: change this so that more errors are detected and so more possible to not
  // return 0
  return 0;
}

template<typename TcalcEngine, typename Tstore, std::size_t numOfRandoms>
void Driver<TcalcEngine, Tstore, numOfRandoms>::SetParameters( Parameters &MyParameters )
{
  // Do specific setup for this driver object
  // The only settings necessary are for the random number generator
  // TODO: clean up the 'Parameters' object, so this doesn't look so convoluted

  std::list<std::string> ReturnedParamsString;
  
  // The type of random number generator
  // Current options are 'Drand48' and 'Sobol'
  ReturnedParamsString.empty();
  ReturnedParamsString = MyParameters.GetParametersString( "@RandomNumberGenerator" );
  
  SwRandoms::RandDrand48 ranDrand48( numOfRandoms );
  SwRandoms::RandSobol ranSobol( numOfRandoms );

  if ( ReturnedParamsString.front() == "Drand48" )
    randomGenerator = ranDrand48;
  else if ( ReturnedParamsString.front() == "Sobol" )
    randomGenerator = ranSobol;
  else
  {
    //return -1;
  }

  // Now that we have set up the type of generator, we need a seed for the
  // generator
  // TODO: add this as an option to the Parameters object, including the time(NULL) option
  //long seed = time(NULL);
  long seed = 2;
  randomGenerator->SetSeed( seed );

  //return 0;
}

template<typename TcalcEngine, typename Tstore, std::size_t numOfRandoms>
void Driver<TcalcEngine, Tstore, numOfRandoms>::RunOneIteration()
{
  // The guts of the Monte Carlo calculation
  // For one iteration, we fetch in a set of random numbers,
  // (these are supplied as a uniform distribution over [0,1)
  // then pass this to the calculation object,
  // (which will transform them to the required limits)
  // then iterate over all the points in the deterministic dimensions,
  // getting the calculation result for each point,
  // and passing this result to the data store.
  // TODO: generalize the iteration over deterministic dimensions to n-dimensions,
  // with n chosen as a template parameter to the driver object

  // First, get the array of random numbers
  // TODO: change the interface to the random number generator
  // to ask / get a fixed size array using boost::array
  // At the moment, the change from a boost::array to a std::vector (MyArray)
  // is cumbersome / a hack
  boost::array<double, numOfRandoms> randomNumbers;
  SwArrays::MyArray randomNumbersTmp( numOfRandoms );
  randomGenerator->GetUniforms( randomNumbersTmp );
  for ( std::size_t i=0; i<numOfRandoms; ++i )
    randomNumbers[i] = randomNumbersTmp[i];

  // Pass the random numbers to the calculator object
  calculator.SetRandoms( randomNumbers );

  // Now iterate over the deterministic dimensions,
  // setting the new coordinate for each point in each dimension,
  // getting the answer and adding it to the data store
  // TODO: improve, to have arbitrary number of dimensions
  std::vector<double> tmp(1);
  for (long dim1=0; dim1<storage.GetLengthDim1(); ++dim1)
  {
    calculator.SetCoord( 1, storage.GetCoordDim1( dim1 ) );

    for (long dim2=0; dim2<storage.GetLengthDim2(); ++dim2)
    {
      calculator.SetCoord( 2, storage.GetCoordDim2( dim2 ) );
      calculator.GetAnswer( tmp );
      storage.AddPoint( dim1, dim2, tmp );
    }
  }
}

template<typename TcalcEngine, typename Tstore, std::size_t numOfRandoms>
void Driver<TcalcEngine, Tstore, numOfRandoms>::SaveResults( std::string outputFile )
{
  // Save the results of our simulation to a permanent file
  // We need to include the state of both the data store and
  // the input runtime parameters

  // First, send the Store2D to file
  // The false here => erase any current contents of the file
  storage.WriteToFile( outputFile, false );

  // Now delegate to the params class to output the settings to the end of the file
  // In this way, we have a record of what the inputs were, and a way to resume and
  // add more statistics if wanted.
  // The true here => append to file
  myParameters.WriteToFile( outputFile, true );
}

#endif
