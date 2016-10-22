/** FactorGraph.cpp --- 
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


#include <cassert>
#include <cstring>
#include "FactorGraph.h"



zzhang::CFactorGraph::CFactorGraph(int NofNodes, int *NofStates)
{
     
     assert(NofNodes > 0);
     m_NofNodes = NofNodes;
     m_NofStates = new int[NofNodes];
     memcpy(m_NofStates, NofStates, sizeof(int) * NofNodes);

     m_bi = new Real*[m_NofNodes];
     m_NodeFactors = std::vector<NodeFactor>(m_NofNodes);
     for(int ni = 0; ni < m_NofNodes; ni++)
     {
	  m_bi[ni] = new Real[m_NofStates[ni]];
	  memset(m_bi[ni], 0, sizeof(Real) * m_NofStates[ni]);
	  m_NodeFactors[ni] = NodeFactor(ni, m_NofStates[ni], m_bi[ni]);
     }
     m_CurrentDecode = new int[NofNodes];
     m_BestDecode = new int[NofNodes];
     BestDecodeV = -ZZHANG_DBL_MAX;
     memset(m_CurrentDecode, 0, sizeof(int) * NofNodes);
     memset(m_BestDecode, 0, sizeof(int) * NofNodes);
     auFactor = NULL;
     srandom(32428328);
     m_verbose = false;
     Evid = std::vector<int> (NofNodes, 0);
     MinDualDecrease = 1e-3;
}


void zzhang::CFactorGraph::AddNodeBelief(int Nid, Real* bi)
{
     assert(Nid < m_NofNodes);
     Real* rbi = m_bi[Nid];
     for(int xi = 0; xi < m_NofStates[Nid]; xi++)
     {
	  (*(rbi++)) += *(bi++);
     }
     m_NodeFactors[Nid].FindLocalMax();
}

std::vector<Real> zzhang::CFactorGraph::GetBelief(int Nid)
{
     assert(Nid < m_NofNodes);
     return std::vector<Real>(m_bi[Nid],
			      m_bi[Nid] + m_NofStates[Nid]);
}

bool zzhang::CFactorGraph::AddSparseEdgeNZ(int ei, int ej, double *data, double *mi, double *mj, int nnz, int *nnzIdx)
{
     assert(ei != ej);
     assert(ei < m_NofNodes && ej < m_NofNodes);
     std::set<int> edge;
     edge.insert(ei);
     edge.insert(ej);
     if(m_FactorId.find(edge) == m_FactorId.end())
     {
	  std::vector<CFactorBase *> Nodes(2);
	  SparseEdgeInternal inEdge;
	  inEdge.ei = ei;
	  inEdge.ej = ej;
	  inEdge.data = data;
	  inEdge.mi = mi;
	  inEdge.mj = mj;
	  inEdge.nnz = nnz;
	  inEdge.nnzIdx = nnzIdx;
	  
	  ExternalData exEdge;
	  exEdge.NofNodes = m_NofNodes;
	  exEdge.NofStates = m_NofStates;
	  Nodes[0] = (CFactorBase* ) &m_NodeFactors[ei];
	  Nodes[1] = (CFactorBase* ) &m_NodeFactors[ej];
	  exEdge.SubFactors = Nodes;
	  SparseEdgeFactor *e = new SparseEdgeNZFactor(&inEdge, &exEdge);

	  m_Factors.push_back(e);
	  return true;
     }
     else{
	  //Not implemented
	  return false;
     }
}

bool zzhang::CFactorGraph::AddSparseEdge(int ei, int ej, double *data, double *mi, double *mj, int nnz, int *nnzIdx)
{
     assert(ei != ej);
     assert(ei < m_NofNodes && ej < m_NofNodes);
     std::set<int> edge;
     edge.insert(ei);
     edge.insert(ej);
     if(m_FactorId.find(edge) == m_FactorId.end())
     {
	  std::vector<CFactorBase *> Nodes(2);
	  SparseEdgeInternal inEdge;
	  inEdge.ei = ei;
	  inEdge.ej = ej;
	  inEdge.data = data;
	  inEdge.mi = mi;
	  inEdge.mj = mj;
	  inEdge.nnz = nnz;
	  inEdge.nnzIdx = nnzIdx;
	  
	  ExternalData exEdge;
	  exEdge.NofNodes = m_NofNodes;
	  exEdge.NofStates = m_NofStates;
	  Nodes[0] = (CFactorBase* ) &m_NodeFactors[ei];
	  Nodes[1] = (CFactorBase* ) &m_NodeFactors[ej];
	  exEdge.SubFactors = Nodes;
	  SparseEdgeFactor *e = new SparseEdgeFactor(&inEdge, &exEdge);

	  m_Factors.push_back(e);
	  return true;
     }
     else{
	  //Not implemented
	  return false;
     }
}


bool zzhang::CFactorGraph::AddEdge(int ei, int ej, double *data)
{
     assert(ei != ej);
     assert(ei < m_NofNodes && ej < m_NofNodes);
     std::set<int> edge;
     edge.insert(ei);
     edge.insert(ej);

     if(m_FactorId.find(edge) == m_FactorId.end())
     {
	  std::vector<CFactorBase *> Nodes(2);
	  EdgeInternal inEdge;
	  inEdge.ei = ei;
	  inEdge.ej = ej;
	  inEdge.data = data;
	  ExternalData exEdge;
	  exEdge.NofNodes = m_NofNodes;
	  exEdge.NofStates = m_NofStates;
	  Nodes[0] = (CFactorBase* ) &m_NodeFactors[ei];
	  Nodes[1] = (CFactorBase* ) &m_NodeFactors[ej];
	  exEdge.SubFactors = Nodes;

	  DenseEdgeFactor* e = new DenseEdgeFactor(&inEdge, &exEdge);

	  m_Factors.push_back(e);
	  return true;
     }
     else{
	  //Not implemented
	  return false;
     }
}



zzhang::FactorGraphDualStore * zzhang::CFactorGraph::StoreDual(){
     FactorGraphDualStore *store = new FactorGraphDualStore(m_NofNodes, static_cast<int>(m_Factors.size()));
     for(int i = 0; i < m_NodeFactors.size(); i++)
     {
	  store->NodeStores[i] = m_NodeFactors[i].Store();
     }
     for(int j = 0; j < m_Factors.size(); j++)
     {
	  store->FactorStores[j] = m_Factors[j]->Store();
     }
     if(auFactor) store->AuctionStore = auFactor->Store();
     store->Evid = Evid;
     return store;
     
}


bool zzhang::CFactorGraph::ReStoreDual(FactorGraphDualStore *store)
{
     if(store == NULL) return false;
     for(int i = 0; i < m_NodeFactors.size(); i++)
     {
	  m_NodeFactors[i].ReStore(store->NodeStores[i]);
     }
     for(int j = 0; j < m_Factors.size(); j++)
     {
	  m_Factors[j]->ReStore(store->FactorStores[j]);
     }
     if(auFactor) auFactor->ReStore(store->AuctionStore);
     Evid = store->Evid;
    return true;
}


void zzhang::CFactorGraph::UpdateMessages()
{
           long start = random();
	   if(m_PrimalDualGap.size() != 0){
		for(int ii = 0; ii < m_PrimalDualGap.size(); ii++)
		{
		     long i = m_PrimalDualGap[ii].first;
		     if((-m_Factors[i]->Primal(m_CurrentDecode)) < 1e-8){
			  continue;
		     }
		     m_Factors[i]->UpdateMessages();
		}
	   }
	   else{
		for(int ii = 0; ii < m_Factors.size(); ii++)
		{
		     long i = (ii + start) % m_Factors.size();
		     m_Factors[i]->UpdateMessages();
		}
	   }
	   if(auFactor)
	   {
		auFactor->Auction();
	   }
	   Dual = 0.0;
	   for(int i = 0; i < m_NofNodes; i++)
	   {
		m_CurrentDecode[i] = m_NodeFactors[i].m_LocalMax;
		Dual += m_bi[i][m_CurrentDecode[i]];
	   }
	   //std::cout << "Dual Here1 " << Dual << std::endl;
	   if(auFactor) Dual += auFactor->SumPrice;
	   
	   //std::cout << "Dual Here2 " << Dual << std::endl;
	   double Primal = Dual;
	   m_PrimalDualGap = std::vector< std::pair<int, double> > (m_Factors.size());
	   for(int i = 0; i <m_Factors.size(); i++)
	   {
		double CPrimal = m_Factors[i]->Primal(m_CurrentDecode);
		Primal += CPrimal;
		m_PrimalDualGap[i] = std::pair<int,double>(i, CPrimal);
	   }
	   std::sort(m_PrimalDualGap.begin(),
		     m_PrimalDualGap.end(),
		     [](const std::pair<int, double>& a,
			const std::pair<int, double>& b){
			  return a.second < b.second;
		     });

	   //std::cout << "Largest Gap: " << m_PrimalDualGap[0].second << std::endl;
	   //std::cout << "Smallest Gap: " << m_PrimalDualGap[m_PrimalDualGap.size() - 1].second << std::endl;
	   
	   
	   if(Primal > BestDecodeV)
	   {
		BestDecodeV = Primal;
		memcpy(m_BestDecode, m_CurrentDecode, sizeof(int) * m_NofNodes);
	   }
	   if(m_verbose) std::cout << "Current Dual " << Dual << " Current Primal " << BestDecodeV ;
}


bool zzhang::CFactorGraph::AddGenericGenericSparseFactor(const std::vector<int> &Nodes, const std::vector<std::vector<int> > &NNZs, double *NNZv){
    std::vector<NodeFactor *> NodeFactors(Nodes.size());
    for(int i = 0; i < Nodes.size(); i++)
    {
        NodeFactors[i] = &m_NodeFactors[Nodes[i]];
    }
    
    zzhang::CGeneralSparseFactor * factor = new zzhang::CGeneralSparseFactor(Nodes,
                                                                             NNZs,
                                                                             NNZv,
                                                                             std::vector<Real *>(0),
                                                                             NodeFactors);
    m_Factors.push_back(factor);
    return true;
    
}
