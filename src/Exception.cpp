// ============================================================================
// Exception.cpp
// ..
// ============================================================================

#include "Exception.h"

ExceptionHandler * ExceptionHandler::istance;

void ExceptionHandler::makeExceptionHandler (RunType runtypeaux) {
  istance = new ExceptionHandler(runtypeaux);
}

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
