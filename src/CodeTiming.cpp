//=============================================================================
//  CodeTiming.cpp
//  ...
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
//=============================================================================


#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include "Exception.h"
#include "CodeTiming.h"
#include "Debug.h"
#include <time.h>
using namespace std;


//=============================================================================
//  CodeTiming::CodeTiming
/// Constructor for CodeTiming class
//=============================================================================
CodeTiming::CodeTiming()
{
  Nblock = 0;
  Nlevelmax = 0;
  ttot = 0.0;
  tstart = clock();
}



//=============================================================================
//  CodeTiming::~CodeTiming
/// CodeTiming destructor
//=============================================================================
CodeTiming::~CodeTiming()
{
}



//=============================================================================
//  CodeTiming::StartTimingSection
/// ..
//=============================================================================
void CodeTiming::StartTimingSection(string newblock, int timing_level)
{
  int iblock;

  // If block not in list, then create new entry to timing block
  if ( blockno.find(newblock) == blockno.end() ) {
    iblock = Nblock;
    blockno[newblock] = iblock;
    block[iblock].timing_level = timing_level;
    block[iblock].block_name = newblock;
    Nblock++;
    Nlevelmax = max(Nlevelmax,timing_level);
  }
  // Else, look-up existing timing block
  else {
    iblock = blockno[newblock];
  }

  // Now record time
  block[iblock].tstart = clock();

  return;
}



//=============================================================================
//  CodeTiming::EndTimingSection
/// ..
//=============================================================================
void CodeTiming::EndTimingSection(string s1)
{
  int iblock;

 // If block not in list, then create new entry to timing block
  if ( blockno.find(s1) == blockno.end() ) {
    cout << "Error : looking for incorrect timing block : " << s1 << endl;
  }
  // Else, look-up existing timing block
  else {
    iblock = blockno[s1];
  }

  block[iblock].tend = clock();
  block[iblock].ttot += (double) (block[iblock].tend - block[iblock].tstart) 
    / (double) CLOCKS_PER_SEC;

  return;
}



//=============================================================================
//  CodeTiming::ComputeTimingStatistics
/// ..
//=============================================================================
void CodeTiming::ComputeTimingStatistics(void)
{
  int iblock;
  int level;
  DOUBLE tcount = 0.0;

  tend = clock();

  ttot = (double) (tend - tstart) / (double) CLOCKS_PER_SEC;

  cout << resetiosflags(ios::adjustfield);
  cout << setiosflags(ios::left);

  cout << "------------------------------------------------------" << endl;
  cout << "Total simulation time : " << ttot << endl;
  cout << "------------------------------------------------------" << endl;

  // Output timing data on each hierarchical level
  for (level=1; level<=Nlevelmax; level++) {
    tcount = 0.0;
    cout << "Level : " << level << endl;
    cout << "Block                         Time        %time" << endl;
    cout << "------------------------------------------------------" << endl;
    for (iblock=0; iblock<Nblock; iblock++) {
      if (block[iblock].timing_level != level) continue;
      tcount += block[iblock].ttot;
      block[iblock].tfraction = block[iblock].ttot / ttot;
      cout << setw(30) << block[iblock].block_name
	   << setw(12) << block[iblock].ttot
	   << setw(12) << 100.0*block[iblock].tfraction << endl;
    }
      cout << setw(30) << "REMAINDER"
	   << setw(12) << ttot - tcount
	   << setw(12) << 100.0*(ttot - tcount)/ttot << endl;
    cout << "------------------------------------------------------" << endl;
  }

  cout << resetiosflags(ios::adjustfield);

  return;
}



//=============================================================================
//  CodeTiming::OutputTimingStatistics
/// ..
//=============================================================================
void CodeTiming::OutputTimingStatistics(void)
{
  return;
}
