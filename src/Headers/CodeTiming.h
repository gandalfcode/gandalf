//=================================================================================================
//  CodeTiming.h
//  Contains class for controlling internal timing routines.
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


#ifndef _CODE_TIMING_H_
#define _CODE_TIMING_H_


#include <map>
#include <string>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include "Precision.h"
using namespace std;



//=================================================================================================
//  Structure TimingBlock
/// \brief  Data structure for holiding code timing information
/// \author D. A. Hubber, G. Rosotti
/// \date   01/10/2013
//=================================================================================================
struct TimingBlock
{
  bool timing_flag;
  int timing_level;
  int Ncalled;
  double tstart_wall;
  double tend_wall;
  double ttot_wall;
  double tfraction_wall;
  DOUBLE ttot;
  DOUBLE tfraction;
  clock_t tstart;
  clock_t tend;
  string block_name;

  TimingBlock()
  {
    timing_flag    = false;
    timing_level   = 0;
    Ncalled        = 0;
    ttot           = 0.0;
    ttot_wall      = 0.0;
    tfraction      = 0.0;
    tfraction_wall = 0.0;
    block_name     = "";
  }

};



//=================================================================================================
//  Class CodeTiming
/// \brief   Internal timing control class
/// \details Internal timing control class
/// \author  D. A. Hubber
/// \date    12/02/2014
//=================================================================================================
class CodeTiming
{
 public:

  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  CodeTiming();
  ~CodeTiming();


  // Other functions
  //-----------------------------------------------------------------------------------------------
  void StartTimingSection(string);
  void EndTimingSection(string);
  void ComputeTimingStatistics(string);
  double WallClockTime(void);


  // CodeTiming class variables
  //-----------------------------------------------------------------------------------------------
  static const int Nblockmax=128;          ///< Max. no. of code timing blocks
  static const int Nlevelmax=8;            ///< Max. no. of timing levels
  int level;                               ///< Current timing level
  int Nlevel;                              ///< No. of timing levels
  //int Nblock;                              ///< No. of timing blocks
  double tstart_wall;                      ///< Start of wall clock timing
  double tend_wall;                        ///< End of wall clock timing
  DOUBLE ttot;                             ///< Total time
  DOUBLE ttot_wall;                        ///< Total wall clock time
  clock_t tstart;                          ///< Start of integer clock
  clock_t tend;                            ///< End of integer clock

  int Nblock[Nlevelmax];                   ///< ..
  int timingOrder[Nlevelmax][Nblockmax];   ///< Order (longest to shortest) of timing blocks
  map<string,int> blockmap[Nlevelmax];     ///< Map of timing block names
  TimingBlock block[Nlevelmax][Nblockmax]; ///< Array of timing blocks
  TimingBlock *levelstack[Nlevelmax];      ///< Stack of current blocks being timed

};
#endif
