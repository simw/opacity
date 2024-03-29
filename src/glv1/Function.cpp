#include "Function.h"

Function::Function( unsigned long Dimensionality )
  : dim( Dimensionality ), xmins( Dimensionality ), xmaxs( Dimensionality )
{

}

Function::~Function()
{

}

void Function::SetLimits( MyArray& xmins_, MyArray& xmaxs_ )
{
  xmins = xmins_;
  xmaxs = xmaxs_;
}

///////////////////////////////////////////////////////////////////

YukawaFunction::YukawaFunction( unsigned long Dimensionality )
  : Function( Dimensionality ), mu2s( Dimensionality )
{

}

YukawaFunction::~YukawaFunction()
{

}

Function* YukawaFunction::clone() const
{
  return new YukawaFunction( *this );
}

void YukawaFunction::function( MyArray& xs, MyArray& res )
{
  for (unsigned long i=0; i<GetDimensionality(); ++i)
  {
    res[i] = 2.*xs[i]*mu2s[i] / pow( xs[i]*xs[i] + mu2s[i], 2 );
  }
}

void YukawaFunction::NormedIntegralInverse( MyArray& uniforms, MyArray& xs )
{
  for (unsigned long i=0; i<GetDimensionality(); ++i)
  {
    xs[i] = (pow(GetXmin(i),2) + mu2s[i]) * (pow(GetXmax(i),2) + mu2s[i]);
    xs[i] /= mu2s[i] + pow(GetXmax(i),2) - (pow(GetXmax(i),2)-pow(GetXmin(i),2))*uniforms[i];
    xs[i] -= mu2s[i];
    xs[i] = sqrt(xs[i]);
  }
}

void YukawaFunction::SetParams( MyArray& mu2s_ )
{
  mu2s = mu2s_;
}

///////////////////////////////////////////////////////////////////

UniformFunction::UniformFunction( unsigned long Dimensionality )
  : Function( Dimensionality )
{

}

UniformFunction::~UniformFunction()
{

}

Function* UniformFunction::clone() const
{
  return new UniformFunction(*this);
}

void UniformFunction::function( MyArray& xs, MyArray& res )
{
  for (unsigned long i=0; i<GetDimensionality(); ++i)
  {
    res[i] = 1. / GetRange(i);
  }
}

void UniformFunction::NormedIntegralInverse( MyArray& uniforms, MyArray& xs )
{
  for (unsigned long i=0; i<GetDimensionality(); ++i)
  {
    xs[i] = GetXmin(i) + uniforms[i]*GetRange(i);
  }
}

void UniformFunction::SetParams( MyArray& maxlens_ )
{

}


///////////////////////////////////////////////////////////////////

UnitStepFunction::UnitStepFunction( unsigned long Dimensionality )
  : Function( Dimensionality ), steps( Dimensionality )
{

}

UnitStepFunction::~UnitStepFunction()
{

}

Function* UnitStepFunction::clone() const
{
  return new UnitStepFunction(*this);
}

void UnitStepFunction::function( MyArray& xs, MyArray& res )
{
  for (unsigned long i=0; i<GetDimensionality(); ++i)
  {
    xs[i] < steps[i]  ? res[i] = 1. / (steps[i]-GetXmin(i))  : res[i] = 0.;
  }
}

void UnitStepFunction::NormedIntegralInverse( MyArray& uniforms, MyArray& xs )
{
  for (unsigned long i=0; i<GetDimensionality(); ++i)
  {
    xs[i] = GetXmin(i) + uniforms[i]*steps[i];
  }
}

void UnitStepFunction::SetParams( MyArray& steps_ )
{
  steps = steps_;
}


///////////////////////////////////////////////////////////////////

ExpDecayFunction::ExpDecayFunction( unsigned long Dimensionality )
  : Function( Dimensionality ), lambdas( Dimensionality )
{

}

ExpDecayFunction::~ExpDecayFunction()
{

}

Function* ExpDecayFunction::clone() const
{
  return new ExpDecayFunction( *this );
}

void ExpDecayFunction::function( MyArray& xs, MyArray& res )
{
  res[0] = 1. / lambdas[0] * exp( -((xs[0]-GetXmin(0)) / lambdas[0]) );
  for (unsigned long i=1; i<GetDimensionality(); ++i)
  {
    res[i] = 1. / lambdas[i] * exp( -(xs[i] - xs[i-1]) / lambdas[i] );
  }
}

void ExpDecayFunction::NormedIntegralInverse( MyArray& uniforms, MyArray& xs )
{
  for (unsigned long i=0; i<GetDimensionality(); ++i)
  {
    xs[i] = (1.-uniforms[i])*exp(-GetXmin(i)/lambdas[i]) + uniforms[i]*exp(-GetXmax(i)/lambdas[i]);
    xs[i] = -lambdas[i] * log( xs[i] );
  }
}

void ExpDecayFunction::SetParams( MyArray& params_ )
{
  lambdas = params_;
}
