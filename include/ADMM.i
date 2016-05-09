%module ADMM

%{
        #define SWIG_FILE_WITH_INIT
        #include "ADMM.hpp"
        #include "CVXcanon/src/Solution.hpp"
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double* INPLACE_ARRAY1, int DIM1) {(double* array, int n)}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* array, int m, int n)}

%include "std_vector.i"
%include "std_map.i"
%include "std_string.i"
%include "std_pair.i"
%include "ADMM.hpp"
%include "typemaps.i"
%include "eigen.i"
%eigen_typemaps(Eigen::MatrixXd)

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
   %template(ArgMap) map<string,Eigen::MatrixXd>;
   %template(ArgdMapVector) vector<map<string,Eigen::MatrixXd> >;
}
%inline%{
        
//helper functions
Eigen::MatrixXd numpyToVector(double *array,int n){
	Eigen::Map<Eigen::MatrixXd> v(array,n,1);
	return v;
}

Eigen::MatrixXd numpyToMatrix(double *array,int m, int n){
	Eigen::Map<Eigen::MatrixXd> mat(array,m,n);
	return mat;
}
%}
