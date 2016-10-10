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


#include <assert.h>
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

  void StartTiming() ;
  void EndTiming() ;

};



//=================================================================================================
//  Class CodeTiming
/// \brief   Internal timing control class
/// \details Internal timing control class
/// \author  D. A. Hubber, G. Rosotti
/// \date    12/02/2014
//=================================================================================================
class CodeTiming
{
#if defined MPI_PARALLEL
  bool CPU_Master;
#endif
 public:

  friend class BlockTimer {
    BlockTimer(string& _block_name, int _level, CodeTiming& _parent, bool delayed_start=false)
     : block_name(_block_name), parent(_parent), level(_level), timing_in_progress(false)
    {
      if (not delayed_start)
        StartTiming() ;
    } ;


    void StartTiming() {
      assert(not timing_in_progress) ;
      parent.StartTimingBlock(level, block_name) ;
      timing_in_progress = true ;
    }

    void EndTiming() {
      parent.EndTimingBlock(level, block_name) ;
      timing_in_progress = false ;
    }

    ~BlockTimer() {
      if (timing_in_progress)
        EndTiming() ;
    } ;

  private:
    string block_name ;
    CodeTiming& parent ;
    int level;
    bool timing_in_progress;
  };

  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  CodeTiming();
  ~CodeTiming();

  // Get a timing object for a new region
  BlockTimer StartNewTimer(string block_name) {
    return BlockTimer(level++, block_name) ;
  }

  BlockTimer NewTimer(string block_name) {
    return BlockTimer(level++, block_name, true);
  }
  void ComputeTimingStatistics(string);


  // Private functions
  //-----------------------------------------------------------------------------------------------
 private:
  void StartTimingBlock(int timing_level, string block_name) {
    __check_levels_are_equal(timing_level, level, block_name) ;
    GetBlock(level++, block_name).StartTiming() ;
  }
  void EndTimingBlock(int timing_level, string block_name) {
    __check_levels_are_equal(timing_level+1, level, block_name) ;
    GetBlock(--level, block_name).EndTiming() ;
  }

  const TimingBlock& GetBlock(int level, string block_name) const {
    return blockmap[level][block_name] ;
  }
  TimingBlock& GetBlock(int level, string block_name) {
    return blockmap[level][block_name] ;
  }

  double WallClockTime(void);

  void __check_levels_are_equal(int, int, string) ;


  // CodeTiming class variables
  //-----------------------------------------------------------------------------------------------
  int level;                               ///< Current timing level
  double tstart_wall;                      ///< Start of wall clock timing
  double tend_wall;                        ///< End of wall clock timing
  DOUBLE ttot;                             ///< Total time
  DOUBLE ttot_wall;                        ///< Total wall clock time
  clock_t tstart;                          ///< Start of integer clock
  clock_t tend;                            ///< End of integer clock

  vector<map<string,TimingBlock> > blockmap;     ///< Map of timing block names
};
#endif
