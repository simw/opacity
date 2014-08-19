#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

/**
	@author Simon Wicks <simon@chiting>
*/
class ProgressBar{
  private:
    double min;
    double max;
    double now;
    int nowPercent;
    
    int width;
    bool update;
  
  public:
    void SetNow(double iNow);
    void PrintPreliminaries();
    void PrintProgress();
    void PrintSpinner();
    void PrintFinal();
    
    ProgressBar(double iMin, double iMax);
    ~ProgressBar();

};

#endif
