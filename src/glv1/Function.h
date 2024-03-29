#ifndef FUNCTION_H
#define FUNCTION_H

#include "../Arrays.h"

using namespace SwArrays;

class Function
{
public:
  Function( unsigned long Dimensionality );
  virtual ~Function();
  virtual Function* clone() const = 0;

  virtual void function( MyArray& xs, MyArray& res ) = 0;
  virtual void NormedIntegralInverse( MyArray& xs, MyArray& res ) = 0;

  virtual void SetLimits( MyArray& xmins_, MyArray& xmaxs_ );
  virtual void SetParams( MyArray& params_ ) = 0;

  inline unsigned long GetDimensionality() const;
  inline double GetXmin( unsigned long i ) const;
  inline double GetXmax( unsigned long i ) const;
  inline double GetRange( unsigned long i ) const;

private:
  unsigned long dim;
  MyArray xmins;
  MyArray xmaxs;

};

unsigned long Function::GetDimensionality() const
{
  return dim;
}

double Function::GetXmin( unsigned long i ) const
{
  return xmins[i];
}

double Function::GetXmax( unsigned long i ) const
{
  return xmaxs[i];
}

double Function::GetRange( unsigned long i ) const
{
  return (xmaxs[i] - xmins[i]);
}

class YukawaFunction : public Function
{
public:
  YukawaFunction( unsigned long Dimensionality );
  virtual ~YukawaFunction();
  virtual Function* clone() const;

  virtual void function( MyArray& xs, MyArray& res );
  virtual void NormedIntegralInverse( MyArray& xs, MyArray& res );

  virtual void SetParams( MyArray& params_);

private:
  MyArray mu2s;
};

class UniformFunction : public Function
{
public:
  UniformFunction( unsigned long Dimensionality );
  virtual ~UniformFunction();
  virtual Function* clone() const;

  virtual void function( MyArray& xs, MyArray& res );
  virtual void NormedIntegralInverse( MyArray& xs, MyArray& res );

  virtual void SetParams( MyArray& params_);

private:

};

class UnitStepFunction : public Function
{
public:
  UnitStepFunction( unsigned long Dimensionality );
  virtual ~UnitStepFunction();
  virtual Function* clone() const;

  virtual void function( MyArray& xs, MyArray& res );
  virtual void NormedIntegralInverse( MyArray& xs, MyArray& res );

  virtual void SetParams( MyArray& params_ );

private:
  MyArray steps;
};

class ExpDecayFunction : public Function
{
public:
  ExpDecayFunction( unsigned long Dimensionality );
  virtual ~ExpDecayFunction();
  virtual Function* clone() const;

  virtual void function( MyArray& xs, MyArray& res );
  virtual void NormedIntegralInverse( MyArray& xs, MyArray& res );

  virtual void SetParams( MyArray& params_);

private:
  MyArray lambdas;
};

#endif
