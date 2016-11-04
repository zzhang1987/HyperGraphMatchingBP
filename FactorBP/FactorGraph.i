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


%typemap(out) std::vector<int>{
     int length = $1.size();
     $result = PyArray_FromDims(1, &length, NPY_INT);
     memcpy(PyArray_DATA((PyArrayObject *) $result),&((*(&$1))[0]),sizeof(int)*length);
 }


//%typemap(in) (int *decode) {
//     if (PyList_Check($input)) {
//	  int size = PyList_Size($input);
//	  Py_ssize_t i = 0;
//	  $1 = (int *)malloc((size) * sizeof(int));
//	  for(int i = 0; i < size; i++)
//	  {
//	       PyObject *s = PyList_GetItem($input,i);
//	       if (!PyInt_Check(s)) {
//		    free($1);
//		    PyErr_SetString(PyExc_ValueError, "List items must be integers");
//		    return NULL;
//	       }
//	       $1[i] = PyInt_AsLong(s);
//	  }
//   }
//   else{
//	  PyErr_SetString(PyExc_ValueError, "Expecting a list");
//	  return NULL;
//   }
     
//}

%typemap(freearg) (int *decode){
     if($1) free($1);
}


%include "BaBTypes.h"
%include "FactorStore.h"
%include "FactorGraphStore.h"
%include "FactorGraph.h"
%include "Factor.h"




