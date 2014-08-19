#ifndef WRAPPER_H
#define WRAPPER_H

/**
  @author Simon Wicks copying out from Joshi 'C++ Design Patterns and Derivatives Pricing'
  A wrapper class to enable easy construction and assignment of new objects, but without
  worrying about deleting them when they go out of scope
  Any classes to be in the 'wrapper' must implement the 'clone' function
**/
template< class T>
class Wrapper
{
public:
  /// Constructor with no parameters - start the pointer as null and zero.
  Wrapper() { DataPtr = 0; }
  
  /// Copy constructor - create a copy of the other object, using the 'clone' function
  Wrapper( const T& inner )
  {
    DataPtr = inner.clone();
  }
  
  /// Destructor - if the pointer is non-null, then delete the object
  ~Wrapper()
  {
    if (DataPtr != 0)
      delete DataPtr;
  }
  
  /// Constructor supplying a 'wrapped' object, clone the inner object
  Wrapper( const Wrapper<T>& original)
  {
    if (original.DataPtr != 0)
      DataPtr = original.DataPtr->clone();
    else
      DataPtr = 0;
  }
  
  /// 
  T* GetPointer()
  { return DataPtr; }

  const T* const GetConstPointer() const
  { return DataPtr; }
  
  Wrapper& operator=( const Wrapper<T>& original )
  {
    if (this != &original)
    {
      if (DataPtr != 0)
        delete DataPtr;
        
      DataPtr = (original.DataPtr !=0) ? original.DataPtr->clone() : 0;
    }
    
    return *this;
  }
  
  T& operator*()
  {
    return *DataPtr;
  }
  
  const T& operator*() const
  {
    return *DataPtr;
  }
  
  const T* const operator->() const
  {
    return DataPtr;
  }
  
  T* operator->()
  {
    return DataPtr;
  }
  
private:
  T* DataPtr;
  
};

#endif
