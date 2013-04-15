//=============================================================================
// Exception.cpp
// ..
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
