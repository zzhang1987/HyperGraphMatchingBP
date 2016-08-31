%module (package="FactorBP") FactorGraph

%{
     #include "BaBTypes.h"
     #include "FactorGraph.h"
     #include "Factor.h"
     #include "FactorStore.h"
     #include "FactorGraphStore.h"
%}

%include "std_vector.i"
%include "std_unordered_map.i"
%include "carrays.i"

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


%include "BaBTypes.h"
%include "FactorStore.h"
%include "FactorGraphStore.h"
%include "FactorGraph.h"
%include "Factor.h"




