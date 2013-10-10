//=============================================================================
//  Exception.cpp
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
//=============================================================================


#include "Exception.h"

ExceptionHandler * ExceptionHandler::istance;


//=============================================================================
//  ExceptionHandler::makeExceptionHandler
/// ...
//=============================================================================
void ExceptionHandler::makeExceptionHandler (RunType runtypeaux) {
  istance = new ExceptionHandler(runtypeaux);
}



//=============================================================================
//  ExceptionHandler::raise
/// ...
//=============================================================================
void ExceptionHandler::raise(string msg) {
  switch (runtype){
  case cplusplus:
    cout << msg << endl;
    exit(-1);
    break;
  case python:
    throw (SerenError (msg));
    break;
  }
}
