/* File : SphSimulation.i */
%module SphSim
%include "std_string.i"
%include "std_map.i"
%include exception.i


%{
#define SWIG_FILE_WITH_INIT
#include "Exception.h"
#include "SphSimulation.h"
#include "Parameters.h"
#include "SimUnits.h"
#include "Sph.h"
#include "SphSnapshot.h"
#include <signal.h>
#include <string>
#include "Precision.h"

void catch_alarm (int SIG) {
signal(SIGINT, catch_alarm);
throw SerenError("CTRL-C received");
}
%}

%exception SphSimulation::Run {
	signal(SIGINT, catch_alarm);
    try{
        $action
    }
    //catch (int e){
    	//printf("Got %i \n", e);
    	//SWIG_exception(SWIG_RuntimeError, "Error error!");
    	//PyErr_SetString(PyExc_KeyboardInterrupt,"You pressed CTRL-C");
    	//return NULL;
    	//exit(0);
    	//SWIG_exception(KeyboardInterrupt, "You pressed CTRL-C");
    //}
    catch (SerenError e){
    	PyErr_SetString(PyExc_KeyboardInterrupt,e.msg.c_str());
    	return NULL;
    }
}

%exception SphSimulation::Setup {
	try{
		$action
	}
	catch (SerenError &e) {
		PyErr_SetString(PyExc_Exception,e.msg.c_str());
		return NULL;
	}
}

%exception Parameters::ReadParamsFile {
	try{
		$action
	}
	catch (SerenError &e) {
		PyErr_SetString(PyExc_Exception,e.msg.c_str());
		return NULL;		
	}
}

%include "numpy.i"
%init %{
import_array();
ExceptionHandler::makeExceptionHandler(python);
%}
%numpy_typemaps(float, NPY_FLOAT, int)
 /* %include <boost_any.i> */
 
 namespace std {
    %template(map_string_int) map<string, int>;
    %template(map_string_string) map<string, string>;
    %template(map_string_float) map<string, float>;
    /*%template(map_string_any) map<string, boost::any>;*/
}

 /* Applies Numpy black magic */
 %apply (float** ARGOUTVIEW_ARRAY1, int *DIM1) {(float** out_array, int* size_array)}

%include "Precision.h"
%include "SphSimulation.h"
%include "Parameters.h"
%include "SimUnits.h"
%include "Sph.h"
%include "SphSnapshot.h"