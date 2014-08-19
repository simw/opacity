//
// Author: Simon Wicks <simonw@phys.columbia.edu>
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <string>

#include "progressbar.h"
#include "timer.h"

#include "driver.h"
#include "glv1/radcalcerwrapper.h"
#include "glv3/glvradiative3.h"
#include "store2d/store.h"

#include "./glv3/qperpcalculator1.h"
#include "./glv3/qperpcalculator3.h"
#include "./glv3/qperpgenerator3.h"

// #define NO_INPUT

#ifdef OPAC1
const std::size_t NUMRANDOMS = 2;
#endif
#ifdef OPAC2
const std::size_t NUMRANDOMS = 4;
#endif
#ifdef OPAC3
const std::size_t NUMRANDOMS = 6;
#endif
#ifdef OPAC4
const std::size_t NUMRANDOMS = 8;
#endif
#ifdef OPAC5
const std::size_t NUMRANDOMS = 10;
#endif
#ifdef OPAC6
const std::size_t NUMRANDOMS = 12;
#endif
#ifdef OPAC7
const std::size_t NUMRANDOMS = 14;
#endif
#ifdef OPAC8
const std::size_t NUMRANDOMS = 16;
#endif
#ifdef OPAC9
const std::size_t NUMRANDOMS = 18;
#endif
#ifdef OPAC10
const std::size_t NUMRANDOMS = 20;
#endif

/** Get input from the user in order to run a Glv calculation
    If we have #define NO_INPUT above then then the user isn't asked,
    it just runs with a default set of parameters
**/
void GetInput( std::string& resume, std::string& inputFile, std::string& outputFile, 
  long& NumberOfIterations );

int main(int argc, char *argv[])
{
  // Sample program doing a Monte Carlo calculation
  // We need to choose the dimensionality of the Monte Carlo at COMPILE time
  // This is in the NUMRANDOMS macros above
  // Note that the calculation engine (eg GLV) may use more than one random number
  // per eg order in opacity, so NUMRANDOMS is a multiple of opacity

  std::cout << std::endl;
  std::cout << "Monte Carlo: Opacity version 3 by Simon Wicks" << std::endl;
  std::cout << "(using " << NUMRANDOMS << " random numbers)." << std::endl << std::endl;

  // The only input needed from the user is:
  // 1) whether we are resuming from a previous run,
  // 2) the filename of the parameters (+data if resuming) to read in
  // 3) the filename to which to write our results
  // 4) the number of iterations to do
  long NumberOfIterations;
  std::string resume = ""; std::string inputFile = ""; std::string outputFile = "";
  GetInput( resume, inputFile, outputFile, NumberOfIterations );

  // The templatized nature of the Driver means that a number of different calculations
  // can be done just by changing the code here in main.
  // A number of examples are given below
  // Just get rid of the /* */ around the block

  /*
  // Glv calculation using most of the old code ('version 1'), but using a wrapper
  // to put it into the new templatized Driver
  // GlvCalcer only really does a dq dphi_q integration and throws away an extra
  // random number per dimension, hence it uses 3 random number per order in opacity
  // ie opacity = NUMRANDOMS/3
  // TODO: change this to work again
  //Driver<RadCalcerWrapper<9>, Store2D, 9> myDriver;
  */

  /*
  // Glv calculation, using the new ('version 3') calculation code (GlvRadiative3),
  // but using a wrapper around the old ('version 1') QperpCalculator code (QperpCalculator1)
  // (which is slower than the newer version, but is a good test to give the same answer)
  // GlvRadiative3 does a dq dphi_q integration, hence opacity = NUMRANDOMS/2
  Driver<GlvRadiative3
    <QperpGenerator3<NUMRANDOMS/2>, QperpCalculator1<NUMRANDOMS/2>, NUMRANDOMS>, 
    Store2D, NUMRANDOMS> myDriver;
  */
  
  // Glv calculation using the new calculation code (GlvRadiative3) and the
  // new qperpcalculator code (QperpCalculator3)
  // GlvRadiative3 does a dq dphi_q integration, hence opacity = NUMRANDOMS/2
  Driver<GlvRadiative3
    <QperpGenerator3<NUMRANDOMS/2>, QperpCalculator3<NUMRANDOMS/2>, NUMRANDOMS>, 
    Store2D, NUMRANDOMS> myDriver;

  /*
  // TODO: write a GlvRadiative4 which includes the z integration
  */


  // End of different calculation examples

  // Call the driver to set up itself and all the other objects
  // If ok, it'll return 0. If not, then it'll return something else.
  int errCode = myDriver.Setup( resume, inputFile );
  if ( errCode != 0 )
  {
    std::cerr << "Error on setting up, code = " << errCode << std::endl;
    return EXIT_SUCCESS;
  }

  // Everything is ready in the calculation, just set up the final
  // bits for the user interface: timer and progress bar
  // ProgressBar: The first number of the current point (start = 0),
  // the second is the end point, the total number
  Timer timer;
  ProgressBar progBar(0,NumberOfIterations);
  progBar.PrintPreliminaries();
  
  // And now we run the main loop
  timer.StartTimer();
  for (long i=0; i<NumberOfIterations; ++i)
  {
    progBar.SetNow(i);
    progBar.PrintProgress();
    myDriver.RunOneIteration();
  }
  progBar.SetNow(NumberOfIterations);
  progBar.PrintProgress();
  timer.StopTimer();
  // And the main loop is done!

  // Tell the driver to save our results to file
  myDriver.SaveResults( outputFile );
 
  // Print to screen a few details of the run
  progBar.PrintFinal();
  double totTime = timer.GetResult();
  std::cout << "Total time = " << totTime << " seconds\n";
  std::cout << "Time per iteration = " << totTime/static_cast<double>(NumberOfIterations);
  std::cout << "\n";
  std::cout << std::endl;

  return EXIT_SUCCESS;
}

void GetInput( std::string& resume, std::string& inputFile, std::string& outputFile, 
  long& NumberOfIterations)
{
  // Get required input for the Monte Carlo run
  // First, ask whether we are resuming the Monte Carlo statistics from a previous run
  do
  {
    std::cout << "Resume from a previous run? (y/n) ";
#ifdef NO_INPUT
   resume = "n";
   std::cout << resume;
#else
    std::cin >> resume;
#endif
    std::cout << "\n";
  } 
  while ( !( resume=="y" || resume=="n" ) );
  
  // Second, filenames
  std::cout << "\nEnter input file and output file\n";
  std::cout << "(no spaces in file names, separate two names with a space): ";
#ifdef NO_INPUT
  //inputFile = "./results/temp.params"; outputFile = "./results/temp.txt";
  inputFile = "./results/brick1.params"; outputFile = "./results/brick1out.txt";
  std::cout << std::endl << inputFile << " " << outputFile;
#else
  std::cin >> inputFile >> outputFile;
#endif
  std::cout << "\n";
  
  // Finally (thirdly), how many Monte Carlo iterations do we want?
  std::cout << "How many Monte Carlo iterations? ";
#ifdef NO_INPUT
  NumberOfIterations = 100000;
  std::cout << NumberOfIterations;
#else
  std::cin >> NumberOfIterations;
#endif
  std::cout << std::endl << std::endl;
}
