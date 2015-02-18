//=============================================================================
//  Debug.h
//  Contains macro definitions of use debug functions.
//  If the appropriate debug level is set in the Makefile, then the 
//  macros are defined to contain code that outputs messages to screen.
//  If not, then functions are undefined and no additional code is included.
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

#if defined(DEBUG1)
#define debug1(x)   cout << x << endl
#else
#define debug1(x)
#endif

#if defined(DEBUG2)
#define debug2(x)   cout << x << endl
#else
#define debug2(x)
#endif
