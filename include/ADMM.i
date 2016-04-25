%module ADMM

%{
        #define SWIG_FILE_WITH_INIT
        #include "ADMM.h"
        #include "CVXcanon/src/Solution.hpp"
        #include "CVXcanon/include/Eigen/Sparse"
        #include "CVXcanon/include/Eigen/Core"
        #include <cstdio>
%}

%include "std_vector.i"
%include "std_map.i"
%include "std_string.i"
%include "std_pair.i"
%include "ADMM.h"


namespace std {
   %template(StringVector) vector< string >;
   %template(StringVector2D) vector<vector< string > >;
   %template(IntVector) vector< int >;
   %template(IntVector2D) vector<vector< int > >;
   %template(IntVector3D) vector<vector<vector< int > > >;
   %template(DoubleVector) vector< double >;
   %template(DoubleVector2D) vector<vector< double > >;
   %template(DoubleVector3D) vector<vector<vector< double > > >;
   %template(LinOpVector) vector< LinOp * >;
   %template(LinOpVector2D) vector<vector< LinOp * > >;
   %template(Solution) map<int,Eigen::MatrixXd>;
   %template(SolutionVector) vector<map<int,Eigen::MatrixXd> >;
   %template(IntPair) std::pair<int,int>;
   %template(PairVector) std::vector<std::pair<int,int> >;
   %template(PairVector2D) std::vector<std::vector<std::pair<int,int> > >;
}

