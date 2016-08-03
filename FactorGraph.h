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
#include <boost/functional/hash.hpp>
#include "PRTypes.h"
#include "Factor.h"
#include "Auction.h"



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
	  /**
	   * Node Factors are stored sepratedly to other factors.
	   */
	  std::vector<NodeFactor> m_NodeFactors;

	  std::unordered_map<std::set<int>, int, boost::hash<std::set<int> > > m_FactorId;

	  /**
	   * Current Decode
	   */
	  int *m_CurrentDecode;
	  /**
	   * Best Decode
	   */
	  int *m_BestDecode;
	  
	  Real BestDecodeV;
	  /**
	   * Special Factor, Auction Factor
	   */

	  CAuctionFactor * auFactor;
	  

     public:

	  void AddAuctionFactor()
	  {
	       auFactor = new CAuctionFactor(m_NofNodes, m_NofStates,
					     m_bi, m_NodeFactors);
	  }
	  
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
	       delete [] m_CurrentDecode;
	       delete [] m_BestDecode;
	       if(auFactor) delete auFactor;
	  }
	  int GetNofNodes(){return m_NofNodes;}

	  /**
	   * Constructor
	   * @param Number of nodes
	   * @param the states of each node
	   */
	  CFactorGraph(int NofNodes, int* NofStates);

	  void AddNodeBelief(int Nid, double* bi);
	  bool AddEdge(int ei, int ej, double *data);
	  bool AddSparseEdge(int ei, int ej, double *data, double *mi, double *mj, int nnz, int *nnzIdx);

	  void UpdateMessages()
	  {
	       for(int i = 0; i < m_Factors.size(); i++)
	       {
		    m_Factors[i]->UpdateMessages();
	       }
	       if(auFactor) auFactor->Auction();
	       double Dual = 0.0;
	       for(int i = 0; i < m_NofNodes; i++)
	       {
		    m_CurrentDecode[i] = m_NodeFactors[i].m_LocalMax;
		    Dual += m_bi[i][m_CurrentDecode[i]];
	       }
	       if(auFactor) Dual += auFactor->SumPrice;
	       double Primal = Dual;
	       for(int i = 0; i <m_Factors.size(); i++)
	       {
		    Primal += m_Factors[i]->Primal(m_CurrentDecode);
	       }
	       if(Primal > BestDecodeV)
	       {
		    BestDecodeV = Primal;
		    memcpy(m_BestDecode, m_CurrentDecode, sizeof(int) * m_NofNodes);
	       }
	       std::cout << "Current Dual " << Dual << " Current Primal " << Primal << std::endl;
	  }
	  
	  /**
	   * For debug in python;
	   */
	  std::vector<Real> GetBelief(int Nid);

	  void PrintFactorInfo(){
	       for(int i = 0; i < m_NofNodes; i++)
	       {
		    m_NodeFactors[i].Print();
	       }
	       for(int i = 0; i < m_Factors.size(); i++)
	       {
		    m_Factors[i]->Print();
	       }
	  }
	  
     };
}

#include "Factor.h"
#endif // FACTORGRAPH_H
