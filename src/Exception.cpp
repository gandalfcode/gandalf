// ============================================================================
// Exception.cpp
// ..
// ============================================================================

#include "Exception.h"

ExceptionHandler * ExceptionHandler::istance;

void ExceptionHandler::makeExceptionHandler (RunType runtypeaux) {
  istance = new ExceptionHandler(runtypeaux);
}
