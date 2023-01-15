/* Daniel R. Reynolds
   SMU Mathematics
   16 September 2015

   This file defines a "fcn" class, that essentially creates a small object to 
   evaluate a user-defined function that carries along its own data parameters. */ 

#ifndef FCN_DEFINED__
#define FCN_DEFINED__


// Fcn class
class Fcn {

public:

  // pure virtual function makes this an 'abstract' base class
  virtual double operator()(double x) = 0;

};

#endif  // FCN_DEFINED__
