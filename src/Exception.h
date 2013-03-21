// ============================================================================
// Exception.h
// ..
// ============================================================================


#ifndef _EXCEPTION_H_
#define _EXCEPTION_H_

#include <cstdlib>
#include <iostream>
#include <string>
using namespace std;


enum RunType {cplusplus, python};

class SerenError {
public:
    SerenError(std::string msgaux){msg = msgaux;}
    ~SerenError() {}
    std::string msg;
};

class StopError: public SerenError {
public:
  StopError(std::string msgaux):
    SerenError(msgaux)
    {};
};

class ExceptionHandler {
private:
  ExceptionHandler(RunType runtypeaux) : runtype(runtypeaux) {}
  static ExceptionHandler * istance;
  const RunType runtype;
public:
  void raise(string msg);
  static void makeExceptionHandler (RunType runtypeaux);
  static ExceptionHandler & getIstance() {return *istance;}
};

#endif // _EXCEPTION_H_
