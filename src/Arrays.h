#ifndef MYARRAYS_H
#define MYARRAYS_H

#include "boost/multi_array.hpp"
#include <valarray>

namespace SwArrays
{

typedef std::vector<double> MyArray;
typedef std::valarray<std::valarray<double> > MyArray2D;

typedef boost::multi_array<double, 1> array1D;
typedef array1D::index index1D;
typedef boost::multi_array<double, 2> array2D;
typedef array2D::index index2D;
typedef boost::multi_array<double, 3> array3D;
typedef array3D::index index3D;

} // End of SwArrays namespace

#endif
