/** FactorGraph.h --- 
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


#ifndef FACTORGRAPH_H
#define FACTORGRAPH_H 1
#include <vector>
#include <unordered_map>
#include "PRTypes.h"
#include "Factor.h"




namespace zzhang{


     
     
     /**
      * A factor graph; This class is used to do MAP inference over an factor graph.
      */
     class CFactorGraph{
    	  
     public:
	  
     private:
	  /**
	   * Number of Nodes;
	   */
	  int m_NofNodes;

	  /**
	   * Number of states for each node.
	   */
	  int *m_NofStates;
	  /**
	   * Higher Factors, 
	   */
	  std::vector< CFactorBase *> m_Factors;
	  /**
	   * Node beliefs
	   */
	  Real**  m_bi;
     public:
	  virtual ~CFactorGraph(){
	       delete [] m_NofStates;
	       for(int i = 0; i < m_Factors.size(); i++)
	       {
		    delete m_Factors[i];
	       }
	       for(int i = 0; i < m_NofNodes; i++)
	       {
		    delete[] m_bi[i];
	       }
	       delete [] m_bi;
	  }
	  int GetNofNodes(){return m_NofNodes;}

	  CFactorGraph(int NofNodes, int* NofStates);

	  void AddNodeBelief(int Nid, Real* bi);

	  /**
	   * For debug in python;
	   */
	  std::vector<Real> GetBelief(int Nid);
	  
     };
}

#include "Factor.h"
#endif // FACTORGRAPH_H
