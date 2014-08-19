
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <string>
#include <boost/array.hpp>

#include "timer.h"

#include "./glv3/qperpcalculator1.h"
#include "./glv3/qperpcalculator3.h"

int main(int argc, char *argv[])
{
  std::cout << "QperpTester" << std::endl;

  boost::array<double, 4> qs; qs[0] = 1; qs[1] = 2; qs[2] = 3; qs[3] = 4;
  boost::array<double, 4> ths; ths[0] = 0; ths[1] = 0; ths[2] = 0; ths[3] = 0;

  QperpCalculator1<4> qp1; QperpCalculator3<4> qp2;
  qp1.SetQsThetas( qs, ths ); qp2.SetQsThetas( qs, ths );
  qp1.SetK( 1. ); qp2.SetK( 1. );

  for (long z=0; z<=7; ++z)
  {
    std::cout << "Qskk[" << z << "]" << std::endl;
    for (long i=1; i<=5; ++i)
    {
      if ( qp1.GetSumQskk( i, z ) != qp2.GetSumQskk( i, z ) )
      {
        std::cout << "(i,z)=(" << i << "," << z << ") ";
        std::cout << "1: " << qp1.GetSumQskk( i, z ) << ", 2: " << qp2.GetSumQskk( i, z ) << std::endl;
      }
    }
    std::cout << std::endl;

    std::cout << "Qs1k[" << z << "]" << std::endl;
    for (long i=1; i<=5; ++i)
    {
      if ( qp1.GetSumQs1k( i, z ) != qp2.GetSumQs1k( i, z ) )
      {
        std::cout << "(i,z)=(" << i << "," << z << ") ";
        std::cout << "1: " << qp1.GetSumQs1k( i, z ) << ", 2: " << qp2.GetSumQs1k( i, z ) << std::endl;
      }
    }
    std::cout << std::endl << std::endl;
  }


  QperpCalculator1<1> qp_v1_1; QperpCalculator3<1> qp_v2_1;
  QperpCalculator1<10> qp_v1_10; QperpCalculator3<10> qp_v2_10;
  Timer timer;
  unsigned long num;

  boost::array<double, 1> qs1, ths1;
  qs1[0] = 1.; ths1[0] = 0.23;
  num = 500000;

  std::cout << "Version 1 (old): " << num << " at n=1: " << std::endl;
  timer.StartTimer();
  for (unsigned long i=0; i!=num; ++i)
    qp_v1_1.SetQsThetas(qs1, ths1);
  timer.StopTimer(); timer.DisplayResult(1); timer.DisplayResult(num);
  std::cout << std::endl;

  std::cout << "Version 2 (new): " << num << " at n=1: " << std::endl;
  timer.StartTimer();
  for (unsigned long i=0; i!=num; ++i)
    qp_v2_1.SetQsThetas(qs1, ths1);
  timer.StopTimer(); timer.DisplayResult(1); timer.DisplayResult(num);
  std::cout << std::endl;

  boost::array<double, 10> qs10, ths10;
  for (long i=0; i<10; ++i)
  {
    qs10[i] = 2.; ths10[i] = 1.;
  }
  num = 100;

  std::cout << "Version 1 (old): " << num << " at n=10: " << std::endl;
  timer.StartTimer();
  for (unsigned long i=0; i!=num; ++i)
    qp_v1_10.SetQsThetas(qs10, ths10);
  timer.StopTimer(); timer.DisplayResult(1); timer.DisplayResult(num);
  std::cout << std::endl;

  std::cout << "Version 2 (new): " << num << " at n=10: " << std::endl;
  timer.StartTimer();
  for (unsigned long i=0; i!=num; ++i)
    qp_v2_10.SetQsThetas(qs10, ths10);
  timer.StopTimer(); timer.DisplayResult(1); timer.DisplayResult(num);
  std::cout << std::endl;

/*

  long NumberOfIterations = 10000;
  std::string resume("n");
  std::string inputFile = "./temp.params";
  std::string outputFile = "./temp.txt";
  
  Driver<GlvRadiative1<9>, Store2D, 9> myDriver;
  //Driver<RadCalcerWrapper<9>, Store2D, 9> myDriver;
  myDriver.Setup( resume, inputFile );

  Timer timer;

  // Set up a progress bar for feedback
  // The first number of the current point (start = 0)
  // The second is the end point, the total number
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

  myDriver.SaveResults( outputFile );
 
  // Print to screen a few details of the run
  progBar.PrintFinal();
  double totTime = timer.GetResult();
  std::cout << "Total time = " << totTime << " seconds\n";
  std::cout << "Time per iteration = " << totTime/static_cast<double>(NumberOfIterations);
  std::cout << "\n";
  std::cout << std::endl;

*/

  return EXIT_SUCCESS;
}
