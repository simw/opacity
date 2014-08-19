//
// C++ Interface: simplecalc3
//
// Description: Example calculation class to fit into the driver class
//
// Author: Simon Wicks <simonw@phys.columbia.edu>
//
#ifndef SIMPLECALC3_H
#define SIMPLECALC3_H

#include <boost/array.hpp>
#include "parameters.h"

/**
  @author Simon Wicks <simonw@phys.columbia.edu>
  Example class to implement the necessary calculation interface to the
  driver class. Whatever the coordinates supplied in the deterministic dimensions
  it just returns the number of random numbers supplied.

  Implementation of other calculations can be based off this interface.
*/
template<std::size_t numOfRandoms>
class SimpleCalc3
{
private:

public:
  SimpleCalc3();

  void SetParameter( Parameters myParams );
  void SetCoord(long dimension, double value);
  void SetRandoms( boost::array<double, numOfRandoms> randoms );
  void GetAnswer( SwArrays::MyArray& answers );
};

template<std::size_t numOfRandoms>
SimpleCalc3<numOfRandoms>::SimpleCalc3()
{

}

template<std::size_t numOfRandoms>
void SimpleCalc3<numOfRandoms>::SetParameters( Parameters myParams )
{

}

template<std::size_t numOfRandoms>
void SimpleCalc3<numOfRandoms>::SetCoord( long dimension, double value)
{
  
}

template<std::size_t numOfRandoms>
void SimpleCalc3<numOfRandoms>::SetRandoms( boost::array<double, numOfRandoms> randoms )
{

}

template<std::size_t numOfRandoms>
void SimpleCalc3<numOfRandoms>::GetAnswer( SwArrays::MyArray& answers )
{
  answers[0] = static_cast<double>(numOfRandoms);
}

#endif
