#include <iostream>
#include <cstdlib>
#include <ctime>

#include "timer.h"

namespace
{
  const double cCpNs = 1.5;
}

unsigned long long int nanotime_ia32(void)
{
     unsigned long long int val;
    __asm__ __volatile__("rdtsc" : "=A" (val) : );
     return(val);
}

Timer::Timer()
{
  // clock_id = CLOCK_REALTIME;
      
  tstart.tv_sec = 0;
  tstart.tv_nsec = 0;
  tstop.tv_sec = 0;
  tstop.tv_nsec = 0;
  
  ia32Start = 0;
  ia32Stop = 0;
}


Timer::~Timer()
{
  
}

void Timer::StartTimer()
{
  //clock_gettime(clock_id, &tstart);

  clock1 = clock();
  
  //ia32Start = nanotime_ia32();
}

void Timer::StopTimer()
{
  //clock_gettime(clock_id, &tstop);
  //_totDiff = double(tstop.tv_sec-tstart.tv_sec) + double(tstop.tv_nsec-tstart.tv_nsec)/1.e9;
  //_totDiffNanSecs = double(tstop.tv_sec-tstart.tv_sec)*1.e9 + double(tstop.tv_nsec-tstart.tv_nsec);

  clock2 = clock();
  clockdiff = static_cast<double>(clock2-clock1)/CLOCKS_PER_SEC;
  
  //ia32Stop = nanotime_ia32();  
  //ia32NanoSecs = int(double(ia32Stop - ia32Start)/cCpNs);
}

void Timer::DisplayResult( double div )
{
  //std::cerr << "Total time: " << _totDiffNanSecs << std::endl;
  std::cerr << "Total time: " << clockdiff / div << std::endl;
  //std::cerr << "Total time: " << (ia32NanoSecs) << std::endl;
}

double Timer::GetResult( )
{
  return clockdiff;
}


