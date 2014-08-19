//
// C++ Interface: sw_templatemeta
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SW_TEMPLATEMETA_H
#define SW_TEMPLATEMETA_H

// to calculate the power
template <long a, long b>
struct TPower
{
   enum { value = a*TPower<a, b-1>::value };
};

template <long a>
struct TPower<a, 0>
{
   enum { value = 1 };
};

#endif

