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
#include <vector>
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

  int Pack(vector<char>&) const ;
  void Unpack(vector<char>::const_iterator&);

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

  class BlockTimer ;

  //=================================================================================================
  //  Class __BlockTimerProxy
  /// \brief   Simply Proxy object to avoid copy-on-return of BlockTimer object (see below).
  /// \details Can only be created directly by CodeTiming object member functions.
  /// \author  R. A. Booth
  /// \date    10/10/2016
  //=================================================================================================
  class __BlockTimerProxy {
  private:
    __BlockTimerProxy(const string& block_name_, CodeTiming* parent_, bool delayed_start_=false)
  : block_name(block_name_), parent(parent_), delayed_start(delayed_start_)
  { } ;

  private:
    friend class BlockTimer;
    friend class CodeTiming;

    string block_name;
    CodeTiming* parent;
    bool delayed_start;
  };

  //=================================================================================================
  //  Class BlockTimer
  /// \brief   Object orientated way of handling timing blocks.
  /// \details In the standard usage pattern this class starts timing on construction and finishes
  ///          timing on destruction, copying the result back to its parent CodeTiming object. As an
  ///          alternative, it is also possible to start and end the timing manually
  /// \author  R. A. Booth
  /// \date    10/10/2016
  //=================================================================================================
  class BlockTimer {
  public:
    BlockTimer(const __BlockTimerProxy& proxy)
    : block_name(proxy.block_name), parent(proxy.parent), level(-1), timing_in_progress(false)
    {
      if (not proxy.delayed_start)
        StartTiming() ;
    } ;

    void StartTiming() {
      assert(not timing_in_progress) ;
      assert(parent != NULL) ;
      level = parent->StartTimingBlock(block_name) ;
      timing_in_progress = true ;
    }

    void EndTiming() {
      parent->EndTimingBlock(level, block_name) ;
      timing_in_progress = false ;
    }

    BlockTimer& operator=(const __BlockTimerProxy& proxy){
      if (timing_in_progress)
        EndTiming() ;
     return *this = BlockTimer(proxy) ;
    }

    ~BlockTimer() {
      if (timing_in_progress)
        EndTiming() ;
    }

    friend class CodeTiming ;

  private:
      BlockTimer& operator=(const BlockTimer& timer) {
        block_name = timer.block_name ;
        parent = timer.parent ;
        level = timer.level ;
        timing_in_progress = timer.timing_in_progress;

        return *this ;
      }

    string block_name ;
    CodeTiming* parent ;
    unsigned int level;
    bool timing_in_progress;
  };


  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  CodeTiming();
  ~CodeTiming();

  // Create a timer for a new section. StartNewTimer also begins the timing.
  //-----------------------------------------------------------------------------------------------
  __BlockTimerProxy StartNewTimer(string block_name);
  __BlockTimerProxy NewTimer(string block_name);

  void ComputeTimingStatistics(string);
  double RunningTime() const;

  float GetBlockTime(string _blockString) {
    std::map<string,TimingBlock>::iterator it;
    for (unsigned int l=0; l<totals.size(); l++) {
      it = totals[l].find(_blockString);
      if (it != totals[l].end()) {
        return (float) it->second.ttot_wall;
      }
    }
    return 0.0f;
  }

  // Private functions
  //-----------------------------------------------------------------------------------------------
 private:
  int StartTimingBlock(const string& block_name) {
    int level = activeblocks.size() ;
    activeblocks.push_back(block_name) ;

    TimingBlock& blk = GetBlock(level, block_name) ;
    blk.block_name = block_name;
    blk.timing_level = level ;
    blk.StartTiming();

    return level ;
  }
  void EndTimingBlock(unsigned int timing_level, const string& block_name) {
    __check_timing_level(timing_level, block_name) ;
    GetBlock(timing_level, block_name).EndTiming() ;
    activeblocks.pop_back();
  }
  TimingBlock& GetBlock(unsigned int level, const string& block_name) {
    if (blockmap.size() == level)
      blockmap.push_back(map<string,TimingBlock>()) ;
    assert(blockmap.size() > level) ;

    return blockmap[level][block_name] ;
  }

  void __check_timing_level(unsigned int, const string&) const;

#ifdef MPI_PARALLEL
  void pack_timing_blocks_for_MPI_send(std::vector<char>&) const ;
  void unpack_timing_blocks_from_MPI_send(std::vector<char>::const_iterator&,
                                          std::vector<map<string, TimingBlock> >&) const ;

#endif

  // CodeTiming class variables
  //-----------------------------------------------------------------------------------------------
  double tstart_wall;                          ///< Start of wall clock timing
  double tend_wall;                            ///< End of wall clock timing
  double ttot;                                 ///< Total time
  double ttot_wall;                            ///< Total wall clock time
  clock_t tstart;                              ///< Start of integer clock
  clock_t tend;                                ///< End of integer clock

  vector<map<string,TimingBlock> > blockmap;   ///< Map of timing block names
  vector<map<string,TimingBlock> > totals;     ///< Timing map containing total time data
  vector<string> activeblocks;                 ///< Block currently being timed at each level
};
#endif
