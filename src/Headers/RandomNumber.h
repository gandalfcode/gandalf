//=================================================================================================
//  RandomNumber.h
//  ..
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=================================================================================================


#ifndef _RANDOM_NUMBER_H_
#define _RANDOM_NUMBER_H_


#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Constants.h"
#include "Precision.h"
using namespace std;



//=================================================================================================
//  Class RandomNumber
/// \brief   Base class for generating random number sequence
/// \details ...
/// \author  D. A. Hubber
/// \date    24/05/2014
//=================================================================================================
class RandomNumber
{
 public:

  RandomNumber() {};
  virtual ~RandomNumber() {};

  virtual int intrand(void) = 0;
  virtual long int longintrand(void) = 0;
  virtual FLOAT floatrand(void) = 0;
  virtual DOUBLE doublerand(void) = 0;
  virtual FLOAT gaussrand(FLOAT, FLOAT) = 0;
  virtual void PrintRandomNumberRange(void) = 0;

};



//=================================================================================================
//  Class XorshiftRand
/// \brief   'xorshift' random number generator
/// \details 'xorshift' random number generator
///          (See Numerical Recipes 3rd Ed. Chap 3, wikipedia 'Xorshift')
/// \author  D. A. Hubber
/// \date    24/05/2014
//=================================================================================================
class XorshiftRand : public RandomNumber
{
public:
 //private:

  // Selected full-period triple (ID:A1 from Numerical Recipes)
  // and MLCG modulo 2^64 mapping (ID:D3 from Numerical Recipes)
  static const unsigned long int a1 = 21;
  static const unsigned long int a2 = 35;
  static const unsigned long int a3 = 4;
  static const unsigned long int amod = 4768777513237032717;
  static const FLOAT invrandmax;


 //public:

  // Internal variables for algorithm
  unsigned long int x;

  // Constructor and destructor
  XorshiftRand(unsigned long int _seed): RandomNumber(), x(_seed)
    {
      for (int k=0; k<10; k++) xorshiftrand();
    };
  //~XorshiftRand() {};


  inline unsigned long int xorshiftrand(void)
  {
    x ^= x >> a1;
    x ^= x << a2;
    x ^= x >> a3;
    return x*amod;
  }

  inline int intrand(void) {return (int) xorshiftrand();}
  inline long int longintrand(void) {return (long int) xorshiftrand();}
  inline FLOAT floatrand(void) {return invrandmax*(FLOAT) xorshiftrand();}
  inline DOUBLE doublerand(void) {return invrandmax*(DOUBLE) xorshiftrand();}

  inline FLOAT gaussrand(FLOAT mean, FLOAT sigma)
  {
    FLOAT U = 0.0;
    FLOAT V = 0.0;
    while (U == 0.0) {
      U = floatrand();
      V = floatrand();
    };
    return sqrt(-2.0*log(U))*cos(2*pi*V);
  }

  void PrintRandomNumberRange(void)
  {
    cout << "Integer range;         min : 1,    max : " << pow(2,64) << endl;
    cout << "Floating point range;  min : 0.0,  max : 1.0" << endl;
    return;
  }

};

// Declare invrandmax constant here (prevents warnings with some compilers)
//const FLOAT XorshiftRand::invrandmax = 1.0/1.84467440737095e19;




//=================================================================================================
//  Class DefaultSystemRand
/// \brief   Default random number class which calls system generator.
/// \details Default random number class which calls system generator.
///          Only used when no other generator is selected in parameters file.
/// \author  D. A. Hubber
/// \date    24/05/2014
//=================================================================================================
class DefaultSystemRand : public RandomNumber
{
 public:

  // Constructor and destructor
  DefaultSystemRand(int _seed) : RandomNumber() {};
  //~DefaultSystemRand();

  inline int intrand(void) {return (int) rand()%RAND_MAX;}
  inline long int longintrand(void) {return (long int) rand()%RAND_MAX;}
  inline FLOAT floatrand(void) {return (FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;}
  inline DOUBLE doublerand(void) {return (DOUBLE)(rand()%RAND_MAX)/(DOUBLE)RAND_MAX;}

  inline FLOAT gaussrand(FLOAT mean, FLOAT sigma)
  {
    FLOAT U = 0.0;
    FLOAT V = 0.0;
    while (U == 0.0) {
      U = floatrand();
      V = floatrand();
    };
    return sqrt(-2.0*log(U))*cos(2*pi*V);
  }

  void PrintRandomNumberRange(void)
  {
    cout << "Integer range;         min : 1,    max : " << RAND_MAX << endl;
    cout << "Floating point range;  min : 0.0,  max : 1.0" << endl;
    return;
  }

};
#endif
