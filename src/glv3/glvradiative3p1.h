//
// C++ Interface: glvradiative3p1
//
// Description: An object to fit into 'Driver' to calculate radiation according to GLV - including the z integration
//
// Author: Simon Wicks <simonw@phys.columbia.edu>
//
#ifndef GLVRADIATIVE3P1_H
#define GLVRADIATIVE3P1_H

#include <vector>
#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
#include <complex>

#include "../parameters.h"
#include "../sw_templatemeta.h"
#include "../constants.h"

// Anon namespace to keep these accessible only in this file
namespace
{

const double overhbarc = 5.076142131979695431;
const double overpi = 0.3183098861837906715;

}

/**
  @author Simon Wicks <simonw@phys.columbia.edu>
  Radiative energy loss from Gyulassy-Levai-Vitev, to multiple orders in opacity. This version has an 'arbitrary' distribution of scattering
  along the z direction - hence has an extra z integration per order in
  opacity compared to basic glvradiative3

  TODO: complete this, so the z distribution is supplied as a template parameter


  Similar to glvradiative3:
  In order to allow for testing different various different methods of calculating the 
  necessary q_perp vectors and their necessary dot products, this class is templatized
  to allow the dropping in at compile time of different classes to do these tasks.

  Template parameters:
  TqperpGenerate: takes uniform random numbers supplied, converts them into numbers
    distributed from 0 to q_max, with a 'weight' given to correspond to the Gyulassy-Wang
    collision distribution (or other). The produced randoms can be uniformly distributed
    from 0 to q_max with Yukawa weights, or Yukawa weighted from 0 to q_max with a weight
    of 1, or some other choice.
    This object has the interface:
      .SetParameters( temperature, jet energy, debye mu )
      .GetQsThetas( boost::array<double, numOfRandoms>& randoms,
         boost::array<double, n>& Qs, boost::array<double, n>& Thetas )
      .GetQsEventWeight( long i )

  TqperpCalculate: takes the 2D q_perp vectors in polar coords, finds all necessary
    dot products needed for the GLV calculation. This has the interface:
      .SetK( double k )
      .SetQsThetas( boost::array<double, n>& Qs, boost::array<double, n>& Thetas )
      double .GetSumQskk( long k, long zeroes )
      double .GetSumQs1k( long k, long zeroes )

  numOfRandoms: the number of random numbers needed to produce the result
    (note that the order in opacity is not numOfRandoms, it's given by n - defined below) 
*/
template<typename TqperpGenerate, typename TqperpCalculate, std::size_t numOfRandoms>
class GlvRadiative3p1
{
private:
  /// The order in opacity
  static const long n = numOfRandoms/3;

  /// Takes given uniform random numbers [0,1), gives q_perp vectors (polar coords)
  TqperpGenerate qpgen;
  /// Takes q_perp vectors, calculates sum q_i.q_j's as necessary
  TqperpCalculate qpcalc;

  // A large number of specific parameters for the (D)GLV calculation
  double _mass2;
  double _mg2;
  double _en;
  double _overEn;
  double _length;
  double _cr;
  double _alphas;
  double _gluonlambda;
  double _temperature;
  double _mu;

  double _dglvBeta;
  double _inteferenceFrac;
  double _x;
  double _overx;
  double _k;

  long _switchkmax;
  bool _abovekmax;

  bool _diffexclude;

  /**
    As the GLV calculation only strictly applies for x<<1 and collinear, eikonal radiation,
    a rough kinematic cutoff has to be applied to the k integration. The Monte Carlo might not
    know about this when supplying (x,k) pairs, so we need to simply return 0 for these.
  **/
  double GetKmax() const;
  /** An interference term. In the correlated case, this is a bunch of cosines. But in this
    special case when the scatters are 'uncorrelated' and distributed exp(-L/lambda), this
    gives a simpler form.
  **/
  double Interference( long m, long zeroes ) const;
  /// The GLV CdotB
  double CdotB( long m, long zeroes ) const;
  /** Calls interference*cdotb multiple times, summing over all m values from
    mmin to n. The choice of mmin depends on whether we want to include the classical
    diffusion terms (which don't contribute to dN/dx, but do alter the shape of dN/dxdk).  
  **/
  double Summation( long mmin, long zeroes ) const;
  /** The ultimate aim of this object - returning dN/dxdk
    Note that x and k are set in the 'SetCoord' function
  **/
  double GetdNdxdk() const;

public:
  GlvRadiative3p1();
  // The 4 following public functions satisfy the interface needed by the Driver object
  /// Getting all the parameters for GLV from the input file
  void SetParameters( Parameters& inParams );
  /** Setting x and k
    Note that setting x requires no recalculation of the q-vectors; setting k does need
    a little recalculation. Hence, iterating over x for fixed k will be more efficient
    than iterating over k for fixed x.
  **/
  void SetCoord(long dimension, double value);
  /** Pass in random number distributed from [0,1), use them to set Q and Theta as necessary.
    This uses the TqperpGenerator object to convert into [0,qmax), then uses
    TqperpCalculator to do all necessary calculation for this set of Qs and Thetas
  **/
  void SetRandoms( boost::array<double, numOfRandoms>& randomNumbers );
  /** The Driver calls GetAnswer, to fill a vector of results. Here, this just inlines
    to call the private function GetdNdxdk()
  **/
  void GetAnswer( std::vector<double>& answers ) const;
};

template<typename TqperpGenerate, typename TqperpCalculate, std::size_t numOfRandoms>
GlvRadiative3p1<TqperpGenerate, TqperpCalculate, numOfRandoms>::GlvRadiative3p1()
  : _x( 0. ), _k( 0. )
{

}

template<typename TqperpGenerate, typename TqperpCalculate, std::size_t numOfRandoms>
void GlvRadiative3p1<TqperpGenerate, TqperpCalculate, numOfRandoms>::SetParameters( Parameters& inParams )
{
  std::vector<double> ReturnedParamsDouble;
  std::vector<long> ReturnedParamsLong;
  std::list<std::string> ReturnedParamsString;
  std::list<std::string>::iterator it;

  // Check that the opacity order in the input file matches
  // that of our compiled executable. The opacity is not a runtime choice,
  // it's a compile time choice - this just makes sure that it matches, and that
  // the opacity is understood by the user
  ReturnedParamsLong = inParams.GetParametersLong( "@opacityOrder" );
  long opac = ReturnedParamsLong[0];
  if ( opac != n )
  {
    // The opacity doesn't match!!
    // Remember that n is defined as a static long above
    std::cerr << "Radcalcer, run time opacity doesn't match that set at compile time.\n";
    std::cerr << "For possible optimization purposes, opacity is chosen at compile time";
    std::cerr << std::endl;
    exit(0);
  }

  // First, we get the parameters of the medium
  // This is: 1st = mu, 2nd = temperature, 3rd = gluon mass, 4th = gluon mean free path
  // Here, we only need the gluon mass
  // the zDist will need the others
  ReturnedParamsDouble = inParams.GetParametersDouble( "@mediumParams" );
  _mg2 = ReturnedParamsDouble[2]*ReturnedParamsDouble[2];
  _temperature = ReturnedParamsDouble[1];
  _gluonlambda = ReturnedParamsDouble[3];
  _mu = ReturnedParamsDouble[0];

  // Next, get the jet flavour
  ReturnedParamsString = inParams.GetParametersString( "@jetFlavour" );
  std::string jetFlavour = ReturnedParamsString.front();
  if ( jetFlavour == "Gluon" )
  {
    _cr = 3.; _mass2 = _mg2;
  }
  else if ( jetFlavour == "Light" )
  {
    _cr = 4./3.; _mass2 = _mg2 / 2;
  }
  else if ( jetFlavour == "Charm" )
  {
    _cr = 4./3.; _mass2 = 1.2*1.2;
  }
  else if ( jetFlavour == "Bottom" )
  {
    _cr = 4./3.; _mass2 = 4.75*4.75;
  }
  else
  {
    std::cerr << "@jetFlavour not understood" << std::endl;
    exit(0);
  }

  // Now, check whether also specifying a jet mass
  if ( inParams.IsParamInList( "@jetMassDirect" ) )
  {
    ReturnedParamsDouble = inParams.GetParametersDouble( "@jetMassDirect" );
    if ( ReturnedParamsDouble.size() > 0 )
      _mass2 = ReturnedParamsDouble[0]*ReturnedParamsDouble[0];
  }

  // Now, the momentum, to give the jet energy
  ReturnedParamsDouble = inParams.GetParametersDouble( "@jetMomentum" );
  double jetMomentum = ReturnedParamsDouble[0];
  // Now calculate and set energy, and mass
  _en = sqrt( _mass2 + jetMomentum*jetMomentum );
  _overEn = 1. / _en;

  // The jet path length in the medium
  ReturnedParamsDouble = inParams.GetParametersDouble( "@pathLength" );
  _length = ReturnedParamsDouble[0];

  // The setting on the k max
  ReturnedParamsLong = inParams.GetParametersLong( "@limitSet" );
  _switchkmax = ReturnedParamsLong[0];

  // The strong coupling, alpha_s
  ReturnedParamsString = inParams.GetParametersString( "@alpha" );
  it = ReturnedParamsString.begin();
  if ( *it == "fixed" )
  {
    ++it;
    _alphas = boost::lexical_cast<double>( *it );
  }
  else
  {
    std::cerr << "Radcalcer, @alpha not understood";
    std::cerr << std::endl;
    exit(0);
  }

  ReturnedParamsString = inParams.GetParametersString( "@incClassicalDiffusion" );
  if ( ReturnedParamsString.front() == "yes" )
    _diffexclude = false;
  else if ( ReturnedParamsString.front() == "no" )
    _diffexclude = true;
  else
  {
    std::cerr << "Radcalcer, @incClassicalDiffusion not understood as yes or no";
    std::cerr << std::endl;
    exit(0);
  }

}

template<typename TqperpGenerate, typename TqperpCalculate, std::size_t numOfRandoms>
void GlvRadiative3p1<TqperpGenerate, TqperpCalculate, numOfRandoms>::SetCoord( long dimension, double value)
{
  // Dimension 1 = k
  // Dimension 2 = x

  switch (dimension)
  {
    case 1:
      _k = value;
      qpcalc.SetK( _k );
      break;

    case 2:
      _x = value;
      _overx = 1. / _x;
      _dglvBeta = _mg2*(1.-_x) + _mass2*_x*_x;
      _inteferenceFrac = 0.5*_overx*_overEn*overhbarc;
      break;

    default:

      break;
  }

  ( _k > GetKmax() ) ? _abovekmax = true : _abovekmax = false;

}

template<typename TqperpGenerate, typename TqperpCalculate, std::size_t numOfRandoms>
void GlvRadiative3p1<TqperpGenerate, TqperpCalculate, numOfRandoms>::SetRandoms( boost::array<double, numOfRandoms>& randoms )
{
  boost::array<double, n> Qs;
  boost::array<double, n> Thetas;

  qpgen.SetParameters( _temperature, _en, _mu );
  qpgen.GetQsThetas( randoms, Qs, Thetas );
  qpcalc.SetQsThetas( Qs, Thetas );
}

template<typename TqperpGenerate, typename TqperpCalculate, std::size_t numOfRandoms>
void GlvRadiative3p1<TqperpGenerate, TqperpCalculate, numOfRandoms>::GetAnswer( std::vector<double>& answers ) const
{
  answers[0] = GetdNdxdk();
}

template<typename TqperpGenerate, typename TqperpCalculate, std::size_t numOfRandoms>
double GlvRadiative3p1<TqperpGenerate, TqperpCalculate, numOfRandoms>::GetKmax() const
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

template<typename TqperpGenerate, typename TqperpCalculate, std::size_t numOfRandoms>
double GlvRadiative3p1<TqperpGenerate, TqperpCalculate, numOfRandoms>::Interference( long m, long zeroes ) const
{
  const double modLen = _length / static_cast<double>(n+1) *_inteferenceFrac;

  double result1, result2;
  double omknLen;
  
  std::complex<double> den = std::complex<double>(1.,0.);
  double qpkk;

  for ( long k=2; k<=m; ++k ) 
  {
    qpkk = qpcalc.GetSumQskk( k, zeroes );
    omknLen = (qpkk+_dglvBeta) * modLen;
    den *= std::complex<double>( 1., omknLen );
  }
  
  result1 = den.real() / ( pow(den.real(),2) + pow(den.imag(),2) );
  
  qpkk = qpcalc.GetSumQskk( 1, zeroes );
  omknLen = (qpkk+_dglvBeta)* modLen;
  den *= std::complex<double>( 1., omknLen );
  
  result2 = den.real() / ( pow(den.real(),2) + pow(den.imag(),2) );
  
  /*
  double denreal = 1.; double denimag = 0.;
  double qpkk;

  for ( long k=2; k<=m; ++k ) 
  {
    qpkk = qpcalc.GetSumQskk( k, zeroes );
    omknLen = (qpkk+_dglvBeta) * modLen;
    denreal -= denimag*omknLen; denimag += denreal*omknLen;
  }  
  result1 = denreal / ( denreal*denreal + denimag*denimag );
  
  qpkk = qpcalc.GetSumQskk( 1, zeroes );
  omknLen = (qpkk+_dglvBeta)* modLen;
  denreal -= denimag*omknLen; denimag += denreal*omknLen;  
  result2 = denreal / ( denreal*denreal + denimag*denimag );
  */
  return ( result1 - result2 );
}

template<typename TqperpGenerate, typename TqperpCalculate, std::size_t numOfRandoms>
double GlvRadiative3p1<TqperpGenerate, TqperpCalculate, numOfRandoms>::CdotB( long m, long zeroes ) const
{
  const double over_qs11 = 1./( _dglvBeta + qpcalc.GetSumQskk( 1, zeroes ) );
  const double over_qsmp1mp1 = 1./( _dglvBeta + qpcalc.GetSumQskk( m+1, zeroes ) );
  const double over_qsmm = 1./( _dglvBeta + qpcalc.GetSumQskk( m, zeroes ) );

  const double qs1m = qpcalc.GetSumQs1k( m, zeroes );
  const double qs1mp1 = qpcalc.GetSumQs1k( m+1, zeroes );

  return ( qs1mp1*over_qsmp1mp1 - qs1m*over_qsmm ) * over_qs11;
}

template<typename TqperpGenerate, typename TqperpCalculate, std::size_t numOfRandoms>
double GlvRadiative3p1<TqperpGenerate, TqperpCalculate, numOfRandoms>::Summation( long mmin, long zeroes ) const
{
  double result = 0.;
  for ( long m=mmin; m<=n; ++m )
  {
    result += Interference( m, zeroes ) * CdotB( m, zeroes );
  }
  return -2.*result;
}

template<typename TqperpGenerate, typename TqperpCalculate, std::size_t numOfRandoms>
double GlvRadiative3p1<TqperpGenerate, TqperpCalculate, numOfRandoms>::GetdNdxdk() const
{
  // First, check whether k is over kmax
  if ( _abovekmax )
  {
    return 0.;
  }
  
  double result = 0.;
  for ( long i=0; i<TPower<2,n>::value-1; ++i )
  {
    // If we want to only look at the quantum source term, exclude the classical
    // diffusion terms, then we start at n. If we want all of it, then we start at 1.
    unsigned long mmin;
    _diffexclude ? mmin = n : mmin = 1;
    result += Summation( mmin, i ) * qpgen.GetQsEventWeight( i );
  }

  result *= _k * _overx * overpi  * _cr *_alphas * 2.;
  result *= pow( _length / _gluonlambda, n ) * 1. / SwUtils::_factorial( n );

  return result;
}


#endif
