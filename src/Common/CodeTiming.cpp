//=================================================================================================
//  CodeTiming.cpp
//  Contains class functions for recording CPU and wall-clock times for
//  marked sections of the code.
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


#include <cstdio>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include "Exception.h"
#include "CodeTiming.h"
#include "Debug.h"
#include "InlineFuncs.h"
using namespace std;

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif
#ifdef _OPENMP
#include "omp.h"
#endif


//=================================================================================================
//  WallClockTime
/// Returns current world clock time (from set point) in seconds
//=================================================================================================
inline double WallClockTime()
{
  struct timeval tm;
  gettimeofday( &tm, NULL );
  return (double)tm.tv_sec + (double)tm.tv_usec / 1000000.0;
}


void TimingBlock::StartTiming()
{
 assert(not timing_flag) ;

 tstart      = clock();
 tstart_wall = WallClockTime();

 timing_flag = true ;
}


void TimingBlock::EndTiming()
{
  assert(timing_flag) ;

  tend = clock();
  tend_wall = WallClockTime();

  ttot += (double) (tend - tstart) / (double) CLOCKS_PER_SEC;
  ttot_wall += (double) (tend_wall - tstart_wall);

  Ncalled++ ;
  timing_flag = false ;
}


//=================================================================================================
//  TimingBlock::Pack
/// Pack the timing data from the block for an MPI send
//=================================================================================================
int TimingBlock::Pack(vector<char>& buffer) const
{
  // First pack in to a temporary
  vector<char> temp ;
  append_bytes(temp, &timing_level) ;
  append_bytes(temp, &Ncalled) ;
  append_bytes(temp, &ttot) ;
  append_bytes(temp, &ttot_wall);

  temp.insert(temp.end(), block_name.begin(), block_name.end());

  int size = temp.size();
  int size_tot = size + 2*sizeof(int) ;

  // Add Header, packed data and footer
  append_bytes(buffer, &size) ;
  buffer.insert(buffer.end(), temp.begin(), temp.end()) ;
  append_bytes(buffer, &size) ;

  return size_tot ;
}

//=================================================================================================
//  TimingBlock::Pack
/// Unpack the timing data for the block from an MPI send
//=================================================================================================
void TimingBlock::Unpack(vector<char>::const_iterator& buffer) {
  int size1 ;
  unpack_bytes(&size1, buffer) ;
  int count = 0;
  count += unpack_bytes(&timing_level, buffer) ;
  count += unpack_bytes(&Ncalled,      buffer) ;
  count += unpack_bytes(&ttot,         buffer) ;
  count += unpack_bytes(&ttot_wall,    buffer) ;

  int name_len = size1 - count ;
  block_name = string(buffer, buffer+name_len) ;
  buffer += name_len ;

  int size2;
  unpack_bytes(&size2, buffer) ;

  if (size1 != size2) {
    stringstream message ;
    message << "Error in unpacking of TimingBlock data. Block name: " << block_name
       <<  ". Expected size2=" << size1 <<". Got"
        << "size2=" << size2 ;
    ExceptionHandler::getIstance().raise(message.str());
  }

}



double CodeTiming::RunningTime() const {
  return WallClockTime() - tstart_wall;
}


//=================================================================================================
//  CodeTiming::CodeTiming
/// Constructor for CodeTiming class
//=================================================================================================
CodeTiming::CodeTiming()
{
  ttot        = 0.0;
  ttot_wall   = 0.0;
  tstart      = clock();
  tstart_wall = WallClockTime();
#if defined MPI_PARALLEL
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  CPU_Master = rank==0 ? true : false;
#endif
}



//=================================================================================================
//  CodeTiming::~CodeTiming
/// CodeTiming destructor
//=================================================================================================
CodeTiming::~CodeTiming()
{
}


//=================================================================================================
//  CodeTiming::StartNewTimer
/// Create and start a new timing section
//=================================================================================================
CodeTiming::__BlockTimerProxy CodeTiming::StartNewTimer(string block_name)
{
   return __BlockTimerProxy(block_name, this) ;
 }

//=================================================================================================
//  CodeTiming::StartNewTimer
/// Create a new timing section without starting it
//=================================================================================================
CodeTiming::__BlockTimerProxy CodeTiming::NewTimer(string block_name)
{
   return __BlockTimerProxy(block_name, this,  true);
}

//=================================================================================================
//  CodeTiming::StartNewTimer
/// Check that the timing block that we are trying to destroy is sensible.
//=================================================================================================
void CodeTiming::__check_timing_level(unsigned int level, const string& block_name) const
{
  if (level != activeblocks.size()-1) {
    stringstream message;
    message << "Error in timing block. Level of timer does not agree with level of the block "
        "being timed: " <<  block_name
        << ".Levels: " << level << " " << activeblocks.size() << "\n" ;
    ExceptionHandler::getIstance().raise(message.str());
  }
  if (block_name != activeblocks[level]) {
    stringstream message;
    message << "Error in timing block. Trying to end block name " << block_name << ", but the"
        << "active block is " << activeblocks[level] << "\n" ;
    ExceptionHandler::getIstance().raise(message.str());
  }
}

struct BlockWallTimeSorter
{
  BlockWallTimeSorter(map<string, TimingBlock>& blocks_)
  : blocks(blocks_)
  { } ;

  bool operator()(const string& key1,const string& key2) {
    return blocks[key1].ttot_wall > blocks[key2].ttot_wall ;
  }
private:
  map<string,TimingBlock>& blocks;
} ;


//=================================================================================================
//  CodeTiming::ComputeTimingStatistics
/// Compute all timing statistics of the code and record statistics into an
/// external file 'run_id.timing'.
//=================================================================================================
void CodeTiming::ComputeTimingStatistics
 (string run_id)                       ///< [in] String i.d. of current simulation
{
  DOUBLE tcount = 0.0;                 // Total time
  DOUBLE tcount_wall = 0.0;            // Total wall-clock time
  string filename;                     // Filename for writing statistics
  string fileend = "timing";           // Filename ending
  ofstream outfile;                    // Output file stream object

  filename  = run_id + "." + fileend;
  tend      = clock();
  tend_wall = WallClockTime();
  ttot      = (double) (tend - tstart) / (double) CLOCKS_PER_SEC;
  ttot_wall = tend_wall - tstart_wall ;
  totals    = blockmap;

  int NOpenMP = 1;
  int Nmpi = 1;
#ifdef _OPENMP
  NOpenMP = omp_get_max_threads();
#endif

#if defined MPI_PARALLEL
  int rank ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &Nmpi);

  // First Gather the total clock cycles
  if (rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, &ttot,      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &ttot_wall, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  }
  else {
    MPI_Reduce(&ttot,      &ttot,      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ttot_wall, &ttot_wall, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  }


  // Now gather the block data from each thread
  // First gather the amount of data to send on rank=0
  std::vector<char> send_buffer ;
  pack_timing_blocks_for_MPI_send(send_buffer) ;
  int size = send_buffer.size() ;

  std::vector<int> data_size(Nmpi) ;
  MPI_Gather(&size, 1, MPI_INT, &(data_size[0]), 1, MPI_INT, 0, MPI_COMM_WORLD) ;

  // Allocate space for receive and set offsets of data
  std::vector<char> recv_buffer ;
  std::vector<int> offsets(Nmpi+1) ;
  if (rank == 0) {
    offsets[0] = 0;
    for (int i=0; i < Nmpi; i++)
      offsets[i+1] = offsets[i] + data_size[i] ;

    recv_buffer.resize(offsets.back()) ;
  }

  // Gather the data
  MPI_Gatherv(&(send_buffer[0]), size, MPI_BYTE,
              &(recv_buffer[0]), &(data_size[0]), &(offsets[0]), MPI_BYTE,
              0, MPI_COMM_WORLD) ;

  // Now add the data to the totals
  if (rank == 0) {
    std::vector<char>::const_iterator buf = recv_buffer.begin() ;
    for (int i=0; i < Nmpi; i++) {
      std::vector<map<string, TimingBlock> > tmp_map ;
      unpack_timing_blocks_from_MPI_send(buf, tmp_map) ;
      assert((buf - recv_buffer.begin()) == offsets[i+1]) ;

      if (i == 0) continue ;

      for (unsigned int l=0; l <tmp_map.size(); ++l) {

        typedef map<string,TimingBlock>::iterator iterator ;

        for(iterator iter=totals[l].begin(); iter != totals[l].end(); ++iter) {
          TimingBlock& block = totals[l][iter->first] ;
          block.ttot += iter->second.ttot ;
          block.ttot_wall = max(block.ttot_wall, iter->second.ttot_wall) ;
        }
      }
    }
  }
  else
    return ;
#endif
  ttot_wall = (double) (tend_wall - tstart_wall);

  int NThread = Nmpi * NOpenMP ;

  outfile.open(filename.c_str());
  outfile << resetiosflags(ios::adjustfield);
  outfile << setiosflags(ios::left);
  std::string dashes(100,'-');
  outfile << dashes << endl;
  outfile << "Total simulation wall clock time : " << ttot_wall << endl;

  outfile << "Threads: Total=" << NThread ;
#ifdef MPI_PARALLEL
  outfile << ", MPI=" << Nmpi ;
#endif
#ifdef _OPENMP
  outfile << ", OpenMP=" << NOpenMP ;
#endif
  outfile << endl ;

  outfile << dashes << endl;

  // Output timing data on each hierarchical level
  //-----------------------------------------------------------------------------------------------
  for (unsigned int l=0; l<totals.size(); l++) {
    tcount = 0.0;
    tcount_wall = 0.0;
    outfile << "Level : " << l << endl;
    outfile << "Block";
    outfile << std::string( 35, ' ' );
    outfile << "Max Wall time" << std::string(2,' ') << "%time";
    outfile << std::string(10,' ') << "Av. CPU Time" << std::string(3,' ') << "%time"<<endl;
    outfile << dashes << endl;

    // First sort by wall time
    std::vector<string> keys ;
    for(map<string,TimingBlock>::iterator iter=totals[l].begin();
        iter != totals[l].end(); ++iter)
      keys.push_back(iter->first) ;

    std::sort(keys.begin(), keys.end(), BlockWallTimeSorter(totals[l])) ;

    for(vector<string>::iterator iter=keys.begin(); iter != keys.end(); ++iter) {
      TimingBlock& block = totals[l][*iter] ;

      tcount      += block.ttot;
      tcount_wall += block.ttot_wall;
      block.tfraction      = block.ttot / ttot;
      block.tfraction_wall = block.ttot_wall / ttot_wall;
      outfile << setw(40) << block.block_name
              << setw(15) << block.ttot_wall
              << setw(15) << 100.0*block.tfraction_wall;
      outfile << setw(15) << block.ttot / NThread
              << setw(15) << 100.0*block.tfraction
              << endl;
    }

    outfile << setw(40) << "REMAINDER"
            << setw(15) << ttot_wall - tcount_wall
            << setw(15) << 100.0*(ttot_wall - tcount_wall)/ttot_wall;

    outfile << setw(15) << (ttot - tcount) / NThread
            << setw(15) << 100.0*(ttot - tcount)/ttot
            << endl;
    outfile << dashes << endl;
  }
  //-----------------------------------------------------------------------------------------------

  outfile << resetiosflags(ios::adjustfield);
  outfile.close();

  return;
}

#ifdef MPI_PARALLEL

//=================================================================================================
//  CodeTiming::pack_timing_blocks_for_MPI_send
/// Fills a vector<char> array with the minimal block data from this processor
//=================================================================================================
void CodeTiming::pack_timing_blocks_for_MPI_send(std::vector<char>& buffer) const
{

  // First pack the number of levels, then for each level pack:
  //   Level id
  //   Number of blocks
  //   Block data
  //   -1
  // Finally, pack a -1.

  buffer.clear() ;
  int Nlevels = blockmap.size() ;
  append_bytes(buffer, &Nlevels) ;


  for (int l=0; l < Nlevels; ++l) {
    append_bytes(buffer, &l) ;

    int Nblocks = blockmap[l].size() ;
    append_bytes(buffer, &Nblocks) ;

    typedef map<string, TimingBlock>::const_iterator iterator;

    for (iterator iter=blockmap[l].begin(); iter != blockmap[l].end(); ++iter)
      iter->second.Pack(buffer) ;

    int minus_one = -1 ;
    append_bytes(buffer, &minus_one) ;
  }

  int minus_one = -1 ;
  append_bytes(buffer, &minus_one) ;

}

//=================================================================================================
//  CodeTiming::unpack_timing_blocks_from_MPI_send
/// Create a new blockmap from the data received in an MPI communication.
//=================================================================================================
void CodeTiming::unpack_timing_blocks_from_MPI_send
(vector<char>::const_iterator& buffer,            ///< [in] Buffer to read data from
vector<map<string, TimingBlock> >& newblockmap)   ///< [out] New block map constructed from data
const
{
  // Unpack the data: First the number of levels, then for each level unpack:
  //   Level id
  //   Number of blocks
  //   Block data
  //   -1
  // Finally, unpack a -1.

  int Nlevels ;
  unpack_bytes(&Nlevels, buffer) ;
  newblockmap.resize(Nlevels) ;

  for (int l=0; l < Nlevels; ++l) {
    int lbuf ;
    unpack_bytes(&lbuf, buffer) ;

    if (lbuf != l) {
      stringstream message;
        message << "Error in unpacking of blockmap. Trying get level " << l << ", but got "
                << lbuf << "\n" ;
        ExceptionHandler::getIstance().raise(message.str());
    }

    int Nblocks;
    unpack_bytes(&Nblocks, buffer) ;

    for (int blk = 0; blk < Nblocks; blk++) {
       TimingBlock block ;
       block.Unpack(buffer) ;
       assert(block.timing_level == l) ;
       newblockmap[l][block.block_name] = block ;
    }

    int minus_one;
    unpack_bytes(&minus_one, buffer) ;

    if (minus_one != -1) {
      string message = "Error in unpacking of blockmap, expected end of level data token";
      ExceptionHandler::getIstance().raise(message);
    }
  }

  int minus_one;
  unpack_bytes(&minus_one, buffer) ;
  if (minus_one != -1) {
    string message = "Error in unpacking of blockmap, expected end of data token";
    ExceptionHandler::getIstance().raise(message);
  }
}
#endif
