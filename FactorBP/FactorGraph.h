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
#include <ctime>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <mutex>
#include <thread>
#include "PRTypes.h"
#include "Factor.h"
#include "FactorGraphStore.h"
#include "Auction.h"
#include "BaBTypes.h"
#include "SubTourFactor.h"

namespace zzhang{
     
     /**
      * A factor graph; This class is used to do MAP inference over an factor graph.
      */
     class CFactorGraph{
    	  
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
      
         std::vector< std::mutex > m_NodeMutexes;
         
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
	  /**
	   * Current Best Primal
	   */
	  Real BestDecodeV;
	  /**
	   * Special Factor, Auction Factor
	   */

	  CAuctionFactor * auFactor;
	  /**
	   * Current Dual Objective
	   */
	  double Dual;

	  /**
	   * Current Evindence
	   */
	  std::vector<int> Evid;

	  std::vector< std::pair<int, double> > m_PrimalDualGap;
	  
     public:

	  void SetVerbose(bool verbose){
	       m_verbose = verbose;
	  }
	  /**
	   * Add an auction factor. 
	   */
	  void AddAuctionFactor()
	  {
	       auFactor = new CAuctionFactor(m_NofNodes, m_NofStates,
					     m_bi, m_NodeFactors, Evid, this);
	  }
	  
	  friend class CAuctionFactor;
	  
	  void AddSubTourFactor(int N, int *Nodes, int *AssignMents){
	       SubTourFactor *subTour = new SubTourFactor(N, Nodes, m_NofStates, AssignMents, m_bi, m_NodeFactors);
	       m_Factors.push_back(subTour);
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
	  bool AddSparseEdgeNZ(int ei, int ej, double *data, double *mi, double *mj, int nnz, int *nnzIdx);
         bool AddGenericGenericSparseFactor(const std::vector<int>& Nodes,
                                            const std::vector< std::vector<int> >& NNZs,
                                            double * NNZv);
	  bool m_verbose;
	  double MinDualDecrease;
	  void SetMinDualDecrease(double EPS){
	       MinDualDecrease = EPS;
	  }
	  void Solve(int MaxIter){
	       m_PrimalDualGap.clear();
	       const clock_t begin_time = clock();
	       double lastDual = 1e20;
	       for(int iter=0; iter < MaxIter; iter++)
	       {
		    if(m_verbose)
			 std::cout << "Iter=" << iter << " ";
		    UpdateMessages();
		    if(m_verbose) std::cout << " Time " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;
		    if(fabs(Dual - BestDecodeV) < 1e-4)
			 break;
		    if(Dual < BestDecodeV)
			 break;
		    if(fabs(Dual - lastDual) < MinDualDecrease)
			 break;
		    lastDual = Dual;
		    
	       }
	  }
         
	  /**
	   * Update messages
	   */
	  void UpdateMessages();
	  
	  /**
	   * For debug in python;
	   */
	  std::vector<Real> GetBelief(int Nid);
	  /**
	   *
	   */
	  double ComputeObj(int* decode){
	       //assert(decode.size() == m_NofNodes);
	       //assert(len == m_NofNodes);
	       Real res = 0.0;
	       for(int i = 0; i < m_NofNodes; i++)
	       {
		    res += m_bi[i][decode[i]];
		    if(auFactor)
		    {
			 res += auFactor->prices[decode[i]];
		    }
	       }

	       for(int i = 0; i < m_Factors.size(); i++)
	       {
		    res += m_Factors[i]->Primal(&decode[0]);
	       }
	       return res;
	  }

	  std::vector<int> GetDecode()
	  {
	       return std::vector<int>(m_BestDecode,
				       m_BestDecode + m_NofNodes);
	  }
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
	  bool SetDecode(int NodeId, int Assigned)
	  {
	       assert(NodeId >= 0);
	       if(Evid[NodeId] > 0 && Assigned > 0)
	       {
		    if(Evid[NodeId] != Assigned) return false;
	       }
	       if(Evid[NodeId] > 0 && Assigned < 0)
	       {
		    if(Evid[NodeId] == -Assigned) return false;
	       }
	       Evid[NodeId] = Assigned;
	       if(Assigned > 0){
		    Assigned--;
		    assert(m_NofStates[NodeId] > Assigned);
		    for(int i = 0; i < m_NofStates[NodeId]; i++)
		    {
			 if(i != Assigned) m_bi[NodeId][i] -= 30;
		    }
	       }
	       else{
		    Assigned = -Assigned;
		    Assigned--;
		    assert(m_NofStates[NodeId] > Assigned);
		    m_bi[NodeId][Assigned] -= 30;
	       }
	       return true;
	  }



	  double DualValue(){return Dual;}
	  double PrimalValue(){return BestDecodeV;}

	  int NofNodes(){return m_NofNodes;}
	  
	  MostFractionalNodes FindMostFracNodes()
	  {
	       MostFractionalNodes res;
	       Real Min_Gap = ZZHANG_DBL_MAX;
	       Real eps = 1e-4;
	       int MaxCnt = -1;
	       for(int i = 0; i < m_NofNodes; i++)
	       {
		    Real MaxValue = m_NodeFactors[i].Dual();
		    Real SecMaxValue = -ZZHANG_DBL_MAX;
		    int MaxVID = m_NodeFactors[i].m_LocalMax;
		    int cMaxCnt = 0;
		    for(int xi = 0; xi < m_NofStates[i]; xi++)
		    {
			 
			 if(xi != MaxVID && m_bi[i][xi] > SecMaxValue)
			 {
			      SecMaxValue = m_bi[i][xi];
			 }
			 if(xi != MaxVID && fabs(m_bi[i][xi] - MaxValue) < eps){
			      cMaxCnt++;
			 }
		    }
		    Real gap = fabs(MaxValue - SecMaxValue);
		    if(cMaxCnt > MaxCnt)
		    {
			 MaxCnt = cMaxCnt;
			 Min_Gap = gap;
			 res.Nodes = i;
			 res.States = MaxVID;
			 res.gap = gap;
		    }
	       }
	       return res;
	  }
	  void ResetMax(){
	       BestDecodeV = -1e20;
	  }
     public:
	  FactorGraphDualStore *StoreDual();
	  bool ReStoreDual(FactorGraphDualStore *store);
     };
}

#include "Factor.h"
#endif // FACTORGRAPH_H
