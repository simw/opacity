#ifndef TIMER_H
#define TIMER_H

#include <ctime>

// Clock cycles per nanosecond (600 MHz)

/**
  @author Simon Wicks <simonw@phys.columbia.edu>
  A simple class to time events
  TODO: improve the implementation, make possible different implementations
  At the moment, there is code for a couple of different ways, but these are
  commented out or not.
*/
class Timer{
  private:
    clockid_t clock_id;
    struct timespec tstart,tstop;
    double _totDiff;
    double _totDiffNanSecs;

    clock_t clock1, clock2;
    double clockdiff;

    unsigned long long int ia32Start;
    unsigned long long int ia32Stop;
    unsigned long long int ia32NanoSecs;
  
  public:
    void StartTimer();
    void StopTimer();
    void DisplayResult( double div = 1. );
    double GetResult();
    
    Timer();
    ~Timer();
};

#endif
