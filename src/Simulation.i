//=============================================================================
//  Simulation.i
//  Contains swig interface file for binding Simulation class to python.
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


%module SphSim
%include <std_string.i>
%include "std_string.i"
%include "std_map.i"
%include "std_list.i"
%include exception.i


%{
#define SWIG_FILE_WITH_INIT
#include <signal.h>
#include <string>
#include "Exception.h"
#include "Simulation.h"
#include "Parameters.h"
#include "RiemannSolver.h"
#include "SimUnits.h"
#include "Sph.h"
#include "SphSnapshot.h"
#include "Precision.h"
#include "Render.h"
#include "SmoothingKernel.h"
#include "UnitInfo.h"
#include "HeaderInfo.h"
#include "CodeTiming.h"

void catch_alarm (int SIG) {
signal(SIGINT, catch_alarm);
throw StopError("CTRL-C received");
}
%}

%exception {
  try {        $action    }
   catch (GandalfError e) {
    	PyErr_SetString(PyExc_Exception,e.msg.c_str());
		return NULL;
   }
}

%exception SimulationBase::InteractiveRun {
	signal(SIGINT, catch_alarm);
	PyThreadState *_save;
    _save = PyEval_SaveThread();
    try{
        $action
        PyEval_RestoreThread(_save);
    }
    catch (StopError e){
    	PyEval_RestoreThread(_save);
    	PyErr_SetString(PyExc_KeyboardInterrupt,e.msg.c_str());
    	return NULL;
    }
    catch (GandalfError e){
    	PyEval_RestoreThread(_save);
    	PyErr_SetString(PyExc_Exception,e.msg.c_str());
		return NULL;
    }
}

%exception SimulationBase::GetParam {
	try{
		$action
	}
	catch (GandalfError e) {
		PyErr_SetString(PyExc_Exception,e.msg.c_str());
		return NULL;
	}
}

%exception SimulationBase::SetParam {
	try{
		$action
	}
	catch (GandalfError e) {
		PyErr_SetString(PyExc_Exception,e.msg.c_str());
		return NULL;
	}
}

%exception SimulationBase::PreSetupForPython {
	try{
		$action
	}
	catch (GandalfError &e) {
		PyErr_SetString(PyExc_Exception,e.msg.c_str());
		return NULL;
	}
}

%exception SimulationBase::SetupSimulation {
	try{
		$action
	}
	catch (GandalfError &e) {
		PyErr_SetString(PyExc_Exception,e.msg.c_str());
		return NULL;
	}
}

%exception SimulationBase::ImportArray {
	try {
		$action
	}
	catch (GandalfError &e) {
		PyErr_SetString(PyExc_Exception,e.msg.c_str());
		return NULL;
	}
}

%exception Parameters::ReadParamsFile {
	try{
		$action
	}
	catch (GandalfError &e) {
		PyErr_SetString(PyExc_Exception,e.msg.c_str());
		return NULL;
	}
}

%exception SphSnapshotBase::ExtractArray {
	try{
		$action
	}
	catch (GandalfError &e) {
		PyErr_SetString(PyExc_Exception,e.msg.c_str());
		return NULL;
	}
}

%exception RenderBase::CreateColumnRenderingGrid {
	signal(SIGINT, catch_alarm);
	PyThreadState *_save;
    _save = PyEval_SaveThread();
    try{
        $action
        PyEval_RestoreThread(_save);
    }
    catch (StopError e){
    	PyEval_RestoreThread(_save);
    	PyErr_SetString(PyExc_KeyboardInterrupt,e.msg.c_str());
    	return NULL;
    }
    catch (GandalfError &e) {
    	PyEval_RestoreThread(_save);
		PyErr_SetString(PyExc_Exception,e.msg.c_str());
		return NULL;
	}    
}

%exception RenderBase::CreateSliceRenderingGrid {
	signal(SIGINT, catch_alarm);
	PyThreadState *_save;
    _save = PyEval_SaveThread();
    try{
        $action
        PyEval_RestoreThread(_save);
    }
    catch (StopError e){
        PyEval_RestoreThread(_save);
    	PyErr_SetString(PyExc_KeyboardInterrupt,e.msg.c_str());
    	return NULL;
    }
    catch (GandalfError &e) {
    	PyEval_RestoreThread(_save);
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
    %template(map_string_float) map<string, double>;
    %template(list_snap_base_pointer) list<SphSnapshotBase* >;
    %template(list_strings) list<string>;
    /*%template(map_string_any) map<string, boost::any>;*/
}

 /* Applies Numpy black magic */
 %apply (float** ARGOUTVIEW_ARRAY1, int *DIM1) {(float** out_array, int* size_array)}
 %apply (float* INPLACE_ARRAY1, int DIM1) {(float* values, int Ngrid)}
 %apply (float* ARGOUT_ARRAY1, int DIM1) {(float* vals, int N)}
 %apply (double* IN_ARRAY1, int DIM1) {(double* input, int size)}

 %apply float& OUTPUT { float& scaling_factor };

 %include "HeaderInfo.h"

%include "Precision.h"
%include "RiemannSolver.h"
%include "Simulation.h"
%include "Parameters.h"
%include "SimUnits.h"
%include "Sph.h"
%include "SphSnapshot.h"
%include "Render.h"
%include "SmoothingKernel.h"
%include "UnitInfo.h"
%include "CodeTiming.h"
