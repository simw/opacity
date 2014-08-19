#include <cmath>
#include <iostream>
#include <complex>
#include <boost/lexical_cast.hpp>

#include "radcalcer.h"
#include "../constants.h"

using namespace SwUtils;

RadCalcer::RadCalcer( Parameters& params_, long opacity_, bool correlated_ )
  : zDist( params_, opacity_ ), qperps( opacity_, correlated_ ), n( opacity_ )
{
  std::vector<double> ReturnedParamsDouble;
  std::vector<long> ReturnedParamsLong;
  std::list<std::string> ReturnedParamsString;
  std::list<std::string>::iterator it;

  // First, we get the parameters of the medium
  // This is: 1st = mu, 2nd = temperature, 3rd = gluon mass, 4th = gluon mean free path
  // Here, we only need the gluon mass
  // the zDist will need the others
  ReturnedParamsDouble = params_.GetParametersDouble( "@mediumParams" );
  _mg = ReturnedParamsDouble[2];

  // Next, get the jet flavour
  ReturnedParamsString = params_.GetParametersString( "@jetFlavour" );
  std::string jetFlavour = ReturnedParamsString.front();
  if ( jetFlavour == "Gluon" )
  {
    _cr = 3.; _mass = _mg;
  }
  else if ( jetFlavour == "Light" )
  {
    _cr = 4./3.; _mass = _mg / sqrt(2);
  }
  else if ( jetFlavour == "Charm" )
  {
    _cr = 4./3.; _mass = 1.2;
  }
  else if ( jetFlavour == "Bottom" )
  {
    _cr = 4./3.; _mass = 4.75;
  }
  else
  {
    std::cerr << "@jetFlavour not understood" << std::endl;
    exit(0);
  }

  // Now, check whether also specifying a jet mass
  ReturnedParamsDouble = params_.GetParametersDouble( "@jetMassDirect" );
  if ( ReturnedParamsDouble.size() > 0 )
    _mass = ReturnedParamsDouble[0];
  
  // Now, the momentum, to give the jet energy
  ReturnedParamsDouble = params_.GetParametersDouble( "@jetMomentum" );
  double jetMomentum = ReturnedParamsDouble[0];
  // Now calculate and set energy, and mass
  _en = sqrt( _mass*_mass + jetMomentum*jetMomentum );

  // The jet path length in the medium
  ReturnedParamsDouble = params_.GetParametersDouble( "@pathLength" );
  _length = ReturnedParamsDouble[0];

  // The setting on the k max
  ReturnedParamsLong = params_.GetParametersLong( "@limitSet" );
  _switchkmax = ReturnedParamsLong[0];

  // The strong coupling, alpha_s
  ReturnedParamsString = params_.GetParametersString( "@alpha" );
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

  ReturnedParamsString = params_.GetParametersString( "@incClassicalDiffusion" );
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


  _correlated = correlated_;
}

RadCalcer::~RadCalcer()
{

}

void RadCalcer::DistributeRandoms( MyArray& Randoms )
{
  MyArray inZs( n );
  MyArray inQs( n ), inThs( n ), temps( n ), mu2s( n );
  MyArray qmins( 0., n ), qmaxs( n ), thmins( 0., n ), thmaxs( 2.*pi, n );

  // First, get an array of uniform random inputs for z positions
  for (long i=0; i<n; ++i)
    inZs[i] = Randoms[i];

  zDist.FindRandomPositions( inZs );
  zDist.GetTempsMu2s( temps, mu2s );
  
  // Now get an array of uniform random inputs for qs,thetas
  // and calculate qmaxs from medium params
  for (long i=0; i<n; ++i)
  {
    inQs[i] = Randoms[i + n];
    inThs[i] = Randoms[i + 2*n];
    qmaxs[i] = sqrt(6.*temps[i]*_en);
  }
  
  qperps.FindRandomQs( inQs, qmins, qmaxs, inThs, thmins, thmaxs, mu2s );
}

void RadCalcer::SetXonly( double x )
{
  _x = x;
  _overx = 1./x;

  _dglv = _mg*_mg*(1.-_x) + _mass*_mass*_x*_x;
  // DeltaZi is in fm, _en is in GeV => need factor hbarc
  _frac = 1. / (2. * _x * _en * hbarc);
}

void RadCalcer::SetKonly( double k)
{
  //k = 1.2;
  qperps.SetK( k );
  _k = k;
}

double RadCalcer::_interference( unsigned long m ) const
{
  double term1 = 0.;
  double term2 = 0.;

  for (unsigned long k=2; k<=m; k++)
  {
    term1 += (qperps.GetSumQiQj( k, k )+_dglv) * zDist.GetDeltaZi( k );
  }
  term2 = term1 + (qperps.GetSumQiQj( 1, 1 )+_dglv) * zDist.GetDeltaZi( 1 );

  double tmp = cos(term1*_frac) - cos(term2*_frac);

  return tmp;
}

double RadCalcer::_interferenceExp( unsigned long m ) const
{
  double result1, result2;
  double omknLen;

  //double realpart = 1.; double imagpart = 0.;

  /*for ( unsigned long k=2; k<=m; ++k )
  {
    omknLen = (qperps.GetSumQiQj( k, k )+_dglv)*_frac * _length / static_cast<double>(n+1);
    realpart = realpart*1. - imagpart*omknLen;
    imagpart = imagpart*1. + realpart*omknLen;
  }

  result1 = realpart / (realpart*realpart + imagpart*imagpart);

  omknLen = (qperps.GetSumQiQj( 1, 1 )+_dglv)*_frac * _length / static_cast<double>(n+1);
  realpart = realpart*1. - imagpart*omknLen;
  imagpart = imagpart*1. + realpart*omknLen;

  result2 = realpart / (realpart*realpart + imagpart*imagpart);*/

  std::complex<double> den = std::complex<double>(1.,0.);
  double qpkk;
  for ( unsigned long k=2; k<=m; ++k ) 
  {
    qpkk = qperps.GetSumQiQj( k, k );
    omknLen = (qpkk+_dglv)*_frac * _length / static_cast<double>(n+1);
    den *= std::complex<double>( 1., omknLen );
  }
  
  result1 = den.real() / ( pow(den.real(),2) + pow(den.imag(),2) );
  
  qpkk = qperps.GetSumQiQj( 1, 1 );
  omknLen = (qpkk+_dglv)*_frac * _length / static_cast<double>(n+1);
  den *= std::complex<double>( 1., omknLen );
  
  result2 = den.real() / ( pow(den.real(),2) + pow(den.imag(),2) );
  
  return ( result1 - result2 );
}

double RadCalcer::_cdotb( unsigned long m ) const 
{
  double qs11 = qperps.GetSumQiQj( 1, 1 );
  double qsmp1mp1 = qperps.GetSumQiQj( m+1, m+1 );
  double qsmm = qperps.GetSumQiQj( m, m );

  double qs1m = qperps.GetSumQiQj( 1, m );
  double qs1mp1 = qperps.GetSumQiQj( 1, m+1 );

  double tmp = ( qs1mp1/(qsmp1mp1+_dglv) - qs1m/(qsmm+_dglv) ) / (qs11+_dglv);
  return tmp;
}

void RadCalcer::GetdNdk2dx( MyArray& results )
{
  // Ok, we're going to calculate dNdxdk, and fill the array 'results' with the answer
  // results[n-1] = the full answer, the others are just intermediate stages

  // What recording scheme are we going by?
  // If we set _maxOpac = n here, it will record the full summation, plus n subsets
  // If we set _maxOpac = 2^n, it will record all possible subsets
  //long _maxOpac = n;
  //long _maxOpac = power(2,n);
  long _maxOpac = 1;

  // If we are out of bounds of the k integral, we just return zero for all the whole result array
  if ( _k >= Getkmax() )
  {
    for (long i=0; i<_maxOpac; ++i)
      results[i] = 0.;
    return;
  }

  // Our intermediate stages, and our final result: result
  double msum, qweight, zweight, coeff, result;
  // Reset to zero, will add on terms to this
  result = 0.;

  // If the medium is correlated, then we need to iterate over all the permutations
  // of the delta functions under the integral.
  // If the medium is uncorrelated, we only need to count the delta functions, not 
  // evaluate separately from each. In this way, we save a whole lot of calculation
  // for the uncorrelated medium.
  // We go from 0<z<2^n to 0<z<n. Oooh, orders of magnitude performance enhancement!  

  if ( _correlated )
  {
    // Iterate over all terms in series from the products with the delta function subtracted off
    // Method: count from 1 up to 2^n, the binary representation of that number - of length n -
    // Gives a unique set of zeroed q's. If the binary digit = 1, then that q is zeroed out
    long imin = 0;
    long imax = power(2,n);
    // opac is to do with what we're putting in the results matrix
    // We will increase it as we go along. Wow, what an explanation.
    // long opac = 1;

    // results[n-1] = 0.;
    for (long i=imax-2; i>=imin; --i)
    {
      // Set the zeroed Qs for this term
      qperps.SetZeroedQs( i );
    
      // Sum over all the c dot b terms
      msum = 0.;
      // What m value do we start at?
      unsigned long mmin;
      // If we want to only look at the quantum source term, exclude the classical
      // diffusion terms, then we start at n. If we want all of it, then we start at 1.
      _diffexclude ? mmin = n : mmin = 1;

      // This sum is defined in GLV II.
      for (long m=mmin; m<=n; ++m)
      {
        // If we are in an uncorrelated medium, then we do the z integrals analytically
        // and sum up the resulting Lorentzians. If not, then we have to do the full sum.

        //msum += -2.*_cdotb( m ) * _interference( m );
        msum += -2.*_cdotb( m ) * _interferenceExp( m);
      }

      // Get the q weight - we have generated random q values, but we allow
      // for generating from a different distribution, and then moving a coefficient
      // into the answer. See reweighted Monte-Carlo integration.
      qweight = qperps.GetQeventWeight();
      // Same for the z weight. Remember that we might not even be evaluating a z integral here.
      zweight = zDist.GetZweight();

      // Now add what we have to what we had from the previous z values
      result += qweight * zweight * msum;

      // If we are at a power of 2, put the result into our result matrix
      //if ( i == imax-power(2,opac) )
      //{
      //  (n-opac)%2 == 0 ? results[opac-1] = result : results[opac-1] = -result;
      //  ++opac;
      //}
    }
    results[0] = result;
    // We have left out a few coeffs that are the same for all evaluations.
    // Now multiply through by them.
    coeff = _k * _overx * overpi  *  _cr*_alphas*2.;
    for (long i=0; i!=_maxOpac; ++i)
    {
      results[i] *= coeff;
    }
    // And we're done!
  }
  else
  {
    // We have a correlated medium. We have n separate contributions (as opposed to the 2^n for
    // the correlated medium). In this case, z will represent how many zeroed q's we have.
    long zmin = 0;
    long zmax = n;
    long combin;

    // Now iterate over the possible zeroed out q's
    for (long z=zmax-1; z>=zmin; --z)
    {
      // The combinatoric factor. How many combinations? = n C z
      combin = _Combinatoric( n, z );
      
      // Set the zeroed Qs for this term
      qperps.SetZeroedQs( z );

      // Sum over all the c dot b terms
      msum = 0.;
      // What m value do we start at?
      unsigned long mmin;
      // If we want to only look at the quantum source term, exclude the classical
      // diffusion terms, then we start at n. If we want all of it, then we start at 1.
      _diffexclude ? mmin = n : mmin = 1;

      // This sum is defined in GLV II.
      for (long m=mmin; m<=n; ++m)
      {
        // If we are in an uncorrelated medium, then we do the z integrals analytically
        // and sum up the resulting Lorentzians. If not, then we have to do the full sum.

        //msum += -2.*_cdotb( m ) * _interference( m );
        msum += -2.*_cdotb( m ) * _interferenceExp( m );
      }

      // Get the q weight - we have generated random q values, but we allow
      // for generating from a different distribution, and then moving a coefficient
      // into the answer. See reweighted Monte-Carlo integration.
      qweight = qperps.GetQeventWeight();
      // Same for the z weight. Remember that we might not even be evaluating a z integral here.
      zweight = zDist.GetZweight();

      // Now add what we have to what we had from the previous z values
      result += qweight * zweight * msum;

      results[n-1-z] = result;
    }
    
    // We have left out a few coeffs that are the same for all evaluations.
    // Now multiply through by them.
    coeff = _k * _overx * overpi  *  _cr*_alphas*2.;
    for (long i=0; i!=_maxOpac; ++i)
    {
      results[i] *= coeff;
    }
    // And we're done!
  }
}

