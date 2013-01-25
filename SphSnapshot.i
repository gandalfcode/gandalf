/* File : SphSnapshot.i */
%module SphSnap

%include "std_string.i"
%include "std_map.i"


%{
#define SWIG_FILE_WITH_INIT
#include "SphSnapshot.h"
%}

%include "numpy.i"
%init %{
import_array();
%}
%numpy_typemaps(float, NPY_FLOAT, int)
 /* %include <boost_any.i> */
 
 namespace std {
    %template(map_string_int) map<string, int>;
    /*%template(map_string_any) map<string, boost::any>;*/
}

 /* Applies Numpy black magic */
 %apply (float** ARGOUTVIEW_ARRAY1, int *DIM1) {(float** out_array, int* size_array)}

%include "SphSnapshot.h"
