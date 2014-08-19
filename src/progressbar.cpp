#include <iostream>

#include "progressbar.h"

using namespace std;

ProgressBar::ProgressBar(double iMin, double iMax)
{
  min = iMin;
  max = iMax;
  now = iMin;
  nowPercent = 0;
  
  width = 45;
  
  update = true;
}


ProgressBar::~ProgressBar()
{
}

void ProgressBar::SetNow(double iNow)
{
  int nowPercentTmp;
  
  if (iNow >= min && iNow <= max)
  {
    now = iNow;
    nowPercentTmp = int((now-min)/(max-min) * 100.);
    if (nowPercentTmp != nowPercent)
    {
      nowPercent = nowPercentTmp;
      update = true;
    }
  }
  else
  {
    std::cerr << "ProgressBar: 'now' is out of range, " << iNow << endl;
  }
}

void ProgressBar::PrintPreliminaries()
{
  std::cerr << "Progress: [";
}

void ProgressBar::PrintProgress()
{
  if (update)
  {  
    int numDash = int(width*nowPercent/100.+0.5);
  
    cerr << "\rProgress: [";
    
    int i = 0;
    for (i=1;i<=numDash;i++)
      cerr << "|";
  
    if (numDash < width)
     for (i=numDash+1; i<width; i++)
        cerr << ".";
  
    cerr << "] ";
    PrintSpinner();
    cerr << " (" << int(nowPercent+0.5) << "%)";
  
    update = false;
  }
}

void ProgressBar::PrintSpinner()
{
  static int num = 1;
  switch (num)
  {
    case 1:
      cerr << "|";
      num++;
      break;
      
    case 2:
      cerr << "/";
      num++;
      break;
      
    case 3:
      cerr << "-";
      num++;
      break;
      
    case 4:
      cerr << "\\";
      num = 1;
      break;
  }
}

void ProgressBar::PrintFinal()
{
  cerr << endl;
}
