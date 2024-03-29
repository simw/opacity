#ifndef RADCALCER_H
#define RADCALCER_H

#include "../Arrays.h"
#include "qperpdist.h"
#include "zcolldist.h"
#include "../parameters.h"

using namespace SwArrays;

class RadCalcer
{
private:
  ZposGenerator zDist;
  QperpGenerator qperps;

  const long n;

  double _x;
  double _overx;
  double _k;

  double _en;
  double _mass;
  double _cr;

  double _temp;
  double _mg;
  double _lambda;
  double _length;

  double _alphas;

  unsigned long _switchkmax;

  double _dglv;
  double _frac;

  bool _correlated;
  bool _diffexclude;

  double _interference( unsigned long m ) const;
  double _interferenceExp( unsigned long m ) const;
  double _cdotb( unsigned long m ) const;

public:
  RadCalcer( Parameters& params_, long opacity_, bool correlated_ );
  ~RadCalcer();

  void DistributeRandoms( MyArray& Randoms );
  //void SetRandomZsandQ();
  void SetXonly( double x );
  void SetKonly( double k );

  inline double Getkmax() const;
  void GetdNdk2dx( MyArray& results );

};

double RadCalcer::Getkmax() const
{
  double kmax;
  switch (_switchkmax)
  {
    case 1:
      // Ivan's first version of k_max
      kmax = _x*_en;
    break;

    case 2:
      // Another version from Ivan's code
      _x <= 0.5 ? kmax = _x*_en : kmax = (1.-_x)*_en;
    break;

    case 3:
      // Magdalena's favourite version of k_max
      kmax = 2.*_x*(1.-_x)*_en;
    break;

    default:
      kmax = -1.;
      break;
  }
  return kmax;
}

#endif
