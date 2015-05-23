// SWIG input file of the 'heat' module

%feature("autodoc","1");

%module(docstring=
"'Native' module",
directors="1",
threads="1"
) Native
%{

#include <string>
#include "Native.h"


%}

// ----------- MODULES UTILISES ------------
%include "std_string.i"

// ----------- NATIVE CLASSES ----------------
%include "Native.h"
