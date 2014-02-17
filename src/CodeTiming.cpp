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
  tstart_wall = getRealTime();
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
  block[iblock].tstart_wall = getRealTime();

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

  block[iblock].tend_wall = getRealTime();
  block[iblock].ttot_wall += 
    (double) (block[iblock].tend_wall - block[iblock].tstart_wall);
  
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
  DOUBLE tcount_wall = 0.0;

  tend = clock();
  tend_wall = getRealTime();

  ttot = (double) (tend - tstart) / (double) CLOCKS_PER_SEC;
  ttot_wall = (double) (tend_wall - tstart_wall);

  cout << resetiosflags(ios::adjustfield);
  cout << setiosflags(ios::left);

  cout << "----------------------------------------------------------------------------" << endl;
  cout << "Total simulation time : " << ttot << "    " << ttot_wall << endl;
  cout << "----------------------------------------------------------------------------" << endl;

  // Output timing data on each hierarchical level
  for (level=1; level<=Nlevelmax; level++) {
    tcount = 0.0;
    tcount_wall = 0.0;
    //cout << "Level : " << level << endl;
    cout << "Block                         Time        %time       Wall time   %time" << endl;
    cout << "----------------------------------------------------------------------------" << endl;
    for (iblock=0; iblock<Nblock; iblock++) {
      if (block[iblock].timing_level != level) continue;
      tcount += block[iblock].ttot;
      tcount_wall += block[iblock].ttot_wall;
      block[iblock].tfraction = block[iblock].ttot / ttot;
      block[iblock].tfraction_wall = block[iblock].ttot_wall / ttot_wall;
      cout << setw(30) << block[iblock].block_name
	   << setw(12) << block[iblock].ttot
	   << setw(12) << 100.0*block[iblock].tfraction
	   << setw(12) << block[iblock].ttot_wall
	   << setw(12) << 100.0*block[iblock].tfraction_wall << endl;
    }
    cout << setw(30) << "REMAINDER"
	 << setw(12) << ttot - tcount
	 << setw(12) << 100.0*(ttot - tcount)/ttot
	 << setw(12) << ttot_wall - tcount_wall
	 << setw(12) << 100.0*(ttot_wall - tcount_wall)/ttot_wall << endl;
    cout << "----------------------------------------------------------------------------" << endl;
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




//=============================================================================
//  CodeTiming::getRealTime
/// ..
//=============================================================================
/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */
#if defined(_WIN32)
#include <Windows.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>	/* POSIX flags */
#include <time.h>	/* clock_gettime(), time() */
#include <sys/time.h>	/* gethrtime(), gettimeofday() */

#if defined(__MACH__) && defined(__APPLE__)
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

#else
#error "Unable to define getRealTime( ) for an unknown OS."
#endif



/**
 * Returns the real time, in seconds, or -1.0 if an error occurred.
 *
 * Time is measured since an arbitrary and OS-dependent start time.
 * The returned real time is only useful for computing an elapsed time
 * between two calls to this function.
 */
double CodeTiming::getRealTime( )
{
#if defined(_WIN32)
  FILETIME tm;
  ULONGLONG t;
#if defined(NTDDI_WIN8) && NTDDI_VERSION >= NTDDI_WIN8
  /* Windows 8, Windows Server 2012 and later. ---------------- */
  GetSystemTimePreciseAsFileTime( &tm );
#else
  /* Windows 2000 and later. ---------------------------------- */
  GetSystemTimeAsFileTime( &tm );
#endif
  t = ((ULONGLONG)tm.dwHighDateTime << 32) | (ULONGLONG)tm.dwLowDateTime;
  return (double)t / 10000000.0;
  
#elif (defined(__hpux) || defined(hpux)) || ((defined(__sun__) || defined(__sun) || defined(sun)) && (defined(__SVR4) || defined(__svr4__)))
  /* HP-UX, Solaris. ------------------------------------------ */
  return (double)gethrtime( ) / 1000000000.0;
  
#elif defined(__MACH__) && defined(__APPLE__)
  /* OSX. ----------------------------------------------------- */
  static double timeConvert = 0.0;
  if ( timeConvert == 0.0 )
    {
      mach_timebase_info_data_t timeBase;
      (void)mach_timebase_info( &timeBase );
      timeConvert = (double)timeBase.numer /
	(double)timeBase.denom /
	1000000000.0;
    }
  return (double)mach_absolute_time( ) * timeConvert;
  
#elif defined(_POSIX_VERSION)
  /* POSIX. --------------------------------------------------- */
#if defined(_POSIX_TIMERS) && (_POSIX_TIMERS > 0)
  {
    struct timespec ts;
#if defined(CLOCK_MONOTONIC_PRECISE)
    /* BSD. --------------------------------------------- */
    const clockid_t id = CLOCK_MONOTONIC_PRECISE;
#elif defined(CLOCK_MONOTONIC_RAW)
    /* Linux. ------------------------------------------- */
    const clockid_t id = CLOCK_MONOTONIC_RAW;
#elif defined(CLOCK_HIGHRES)
    /* Solaris. ----------------------------------------- */
    const clockid_t id = CLOCK_HIGHRES;
#elif defined(CLOCK_MONOTONIC)
    /* AIX, BSD, Linux, POSIX, Solaris. ----------------- */
    const clockid_t id = CLOCK_MONOTONIC;
#elif defined(CLOCK_REALTIME)
    /* AIX, BSD, HP-UX, Linux, POSIX. ------------------- */
    const clockid_t id = CLOCK_REALTIME;
#else
    const clockid_t id = (clockid_t)-1;	/* Unknown. */
#endif /* CLOCK_* */
    if ( id != (clockid_t)-1 && clock_gettime( id, &ts ) != -1 )
      return (double)ts.tv_sec +
	(double)ts.tv_nsec / 1000000000.0;
    /* Fall thru. */
  }
#endif /* _POSIX_TIMERS */
  
  /* AIX, BSD, Cygwin, HP-UX, Linux, OSX, POSIX, Solaris. ----- */
  struct timeval tm;
  gettimeofday( &tm, NULL );
  return (double)tm.tv_sec + (double)tm.tv_usec / 1000000.0;
#else
  return -1.0;		/* Failed. */
#endif
}
