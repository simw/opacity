//
// C++ Interface: radcalcerwrapper
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef RADCALCERWRAPPER_H
#define RADCALCERWRAPPER_H

#include "radcalcer.h"
#include "../Arrays.h"
#include "../parameters.h"

using namespace SwArrays;

/**
	@author 
*/
template<std::size_t numOfRandoms>
class RadCalcerWrapper
{
private:
  RadCalcer* radcalc;

public:
  RadCalcerWrapper()
  {

  };

  ~RadCalcerWrapper()
  {
    delete radcalc;
  };

  void SetParameters( Parameters& inParams )
  {
    long opacity = numOfRandoms / 3;
    radcalc = new RadCalcer( inParams, opacity, true );
  };

  void SetCoord(long dimension, double value)
  {
    // Dimension 1 = k
    // Dimension 2 = x

    switch (dimension)
    {
      case 1:
        radcalc->SetKonly( value );
        break;

      case 2:
        radcalc->SetXonly( value );
        break;

      default:
        break;
    }
  };

  void SetRandoms( SwArrays::MyArray& randomNumbers )
  {
    radcalc->DistributeRandoms( randomNumbers );
  };
  
  void GetAnswer( SwArrays::MyArray& answers )
  {
    radcalc->GetdNdk2dx( answers );
  };

};

#endif
