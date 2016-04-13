%module ADMM

%{
        #define SWIG_FILE_WITH_INIT
        #include "ADMM.h"
%}

%include "CVXcanon/src/python/numpy.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_string.i"
%include "ADMM.h"

%init %{
        import_array();
%}

namespace std {
   %template(LinOpVector) vector< LinOp * >;
   %template(LinOpVector2D) vector<vector< LinOp * > >;
}

