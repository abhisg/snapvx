%module ADMM

%{
        #define SWIG_FILE_WITH_INIT
        #include "ADMM.h"
        #include "CVXcanon/src/Solution.hpp"
        #include "CVXcanon/include/Eigen/Sparse"
        #include "CVXcanon/include/Eigen/Core"
        #include <cstdio>
%}

%include "CVXcanon/src/python/numpy.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_string.i"
%include "std_pair.i"
%include "ADMM.h"


%init %{
        import_array();
%}

namespace std {
   %template(IntVector) vector< int >;
   %template(IntVector2D) vector<vector< int > >;
   %template(IntVector3D) vector<vector<vector< int > > >;
   %template(LinOpVector) vector< LinOp * >;
   %template(LinOpVector2D) vector<vector< LinOp * > >;
   %template(solution) map<int,Eigen::MatrixXd>;
   %template(solutionVector) vector<map<int,Eigen::MatrixXd> >;
   %template(IntPair) std::pair<int,int>;
   %template(PairVector) std::vector<std::pair<int,int> >;
}

