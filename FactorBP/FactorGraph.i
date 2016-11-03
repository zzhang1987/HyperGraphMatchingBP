%module (package="FactorBP") FactorGraph
%include "std_vector.i"
%include "carrays.i"

%{
     #define SWIG_FILE_WITH_INIT
     #include "BaBTypes.h"
     #include "FactorGraph.h"
     #include "Factor.h"
     #include "FactorStore.h"
     #include "FactorGraphStore.h"
%}

%include "numpy.i"

%init %{
     import_array();
%}

%array_class(int, intArray);
%array_class(double, doubleArray);




namespace std{
     /* On a side note, the names VecDouble and VecVecdouble can be changed, but the order of first the inner vector matters !*/
  %template(VecInt) vector<int>;
  %template(VecVecInt) vector< vector<int> >;
  %template(VecDouble) vector<double>;
  %template(VecVecdouble) vector< vector<double> >;
 }
%typemap(newfree) FactorGraphDualStore * "delete $1;";
%newobject zzhang::CFactorGraph::StoreDual;

%apply (double * IN_ARRAY1, int DIM1) {(double* seq, int n)};
%apply (int * IN_ARRAY1) {(int* seq)};
%typemap(out) std::vector<int>{
     int length = $1.size();
     $result = PyArray_FromDims(1, &length, NPY_INT);
     memcpy(PyArray_DATA((PyArrayObject *) $result),&((*(&$1))[0]),sizeof(int)*length);
 }



%include "BaBTypes.h"
%include "FactorStore.h"
%include "FactorGraphStore.h"
%include "FactorGraph.h"
%include "Factor.h"




