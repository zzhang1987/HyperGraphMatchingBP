/** Factor.h --- 
 *
 * Copyright (C) 2016 Zhen Zhang
 *
 * Author: Zhen Zhang <zhen@zzhang.org>
 *
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


#ifndef FACTOR_H
#define FACTOR_H 1

#include <utility>
#include <unordered_map>
#include <vector>
#include "PRTypes.h"

#define FACTOR_INVALID -1
#define FACTOR_EDGE_ID 1
#define FACTOR_SPARSEEDGE_ID 2
#define FACTOR_GENERAL_ID 3
#define FACTOR_GENERALSPARSE_ID 4

struct EdgeCreatorStruct{
     int ei;
     int ej;
     Real* data;
     Real* bi;
     Real* bj;
     int* NofStates;
};

struct SparseEdgeCreatorStruct{
     int ei;
     int ej;
     Real *data;
     Real *mi;
     Real *mj;
     int nnz;
     int *nnzIdx;
     Real* bi;
     Real* bj;
     int *NofStates;
};
namespace zzhang{
     /**
      * Base class for a factor.
      */
     class CFactorBase{
     public:
	  /**
	   * A public interface to create a factor.
	   */
	  typedef CFactorBase* (*FactorCreator)(const std::vector<int>& nodes,
						void *data);

	       
	  /**
	   * Register a factor creator.
	   */
	  static bool RegisterFactorCreator(int ID, FactorCreator creator){
	       if(ID == FACTOR_INVALID) return false;
	       if(FactorCreators.find(ID) != FactorCreators.end())
		    return false;
	       FactorCreators[ID] = creator;
	       return true;
	  }

	  /**
	   * Create a factor 
	   */ 
	  static CFactorBase* CreateFactor(int ID, const std::vector<int>& nodes,
					   void *data){
	       if(ID == FACTOR_INVALID) return NULL;
	       if(FactorCreators.find(ID) == FactorCreators.end())
		    return NULL;
	       FactorCreator creator = FactorCreators[ID];
	       return creator(nodes, data);
	  }
	  
	  
	  /**
	   * return the primal value of current factor with given decode.
	   * @param decode given decode
	   */
	  virtual Real Primal(int *decode) = 0;
	  /**
	   * return current dual.
	   */
	  virtual Real Dual() = 0;
	  /**
	   * Updating all message variables.
	   */
	  virtual void UpdateMessages() = 0;

	  virtual ~CFactorBase() = 0;
     private:
	  static std::unordered_map<int, FactorCreator > FactorCreators;
     public:
	  static const int FactorID = FACTOR_INVALID;
     };

     class DenseEdgeFactor : public CFactorBase
     {
     private:
	  Real* bi;
	  Real* bj;
	  Real* bij;
	  int ei;
	  int ej;
	  int *NofStates;
     public:
	  DenseEdgeFactor(EdgeCreatorStruct& edges);
	  virtual Real Primal(int *decode);
	  virtual Real Dual();
	  static const int FactorID = FACTOR_EDGE_ID;
     };
}

#endif // FACTOR_H
