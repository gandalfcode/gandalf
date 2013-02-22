#include <cstdlib>
#include <iostream>
#include <string>
using namespace std;

#ifndef _EXCEPTION_H_
#define _EXCEPTION_H_


enum RunType {cplusplus, python};

class SerenError {
public:
    SerenError(std::string msgaux){msg = msgaux;}
    ~SerenError() {}
    std::string msg;
};

class ExceptionHandler {
private:
  ExceptionHandler(RunType runtypeaux) : runtype(runtypeaux) {}
  static ExceptionHandler * istance;
  const RunType runtype;
public:
  void raise(string msg) {
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
  static void makeExceptionHandler (RunType runtypeaux);
  static ExceptionHandler & getIstance() {return *istance;}
};

#endif // _EXCEPTION_H_
