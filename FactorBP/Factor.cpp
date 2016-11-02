/** Factor.cpp --- 
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
#include "Factor.h"

std::unordered_map<int, zzhang::CFactorBase::FactorCreator> zzhang::CFactorBase::FactorCreators =
     std::unordered_map<int, zzhang::CFactorBase::FactorCreator>();


zzhang::DenseEdgeFactor::DenseEdgeFactor(const void* InParam, const ExternalData* OuParam)
{
     
     assert(OuParam->SubFactors.size() == 2);
     EdgeInternal *internal = (EdgeInternal *) InParam;
     NofStates = OuParam->NofStates;
     n1 = (NodeFactor *)OuParam->SubFactors[0];
     n2 = (NodeFactor *)OuParam->SubFactors[1];
     bi = n1->m_bi; bj = n2->m_bi;
     
     ei = internal->ei; ej = internal->ej;
     int xijMax = NofStates[ei] * NofStates[ej];
     bij = new Real[xijMax];
     memcpy(bij, internal->data, sizeof(Real) * xijMax);
}

zzhang::FactorStore* zzhang::NodeFactor::Store(){
     zzhang::FactorStore *store = new zzhang::FactorStore(m_NofStates);
     memcpy(store->data, m_bi, sizeof(Real) * m_NofStates);
     return store;
}

bool zzhang::NodeFactor::ReStore(zzhang::FactorStore *store)
{
     if(!store) return false;
     memcpy(m_bi, store->data, sizeof(Real) * m_NofStates);
     return store;
}

zzhang::FactorStore* zzhang::SparseEdgeFactor::Store(){
     zzhang::FactorStore *store = new zzhang::FactorStore(NofStates[ei] + NofStates[ej]);
     memcpy(store->data, mi, sizeof(Real) * NofStates[ei]);
     memcpy(store->data + NofStates[ei], mj, sizeof(Real) * NofStates[ej]);
     return store;
}

bool zzhang::SparseEdgeFactor::ReStore(zzhang::FactorStore *store)
{
     if(!store) return false;
     memcpy(mi, store->data, sizeof(Real) * NofStates[ei]);
     memcpy(mj, store->data + NofStates[ei], sizeof(Real) * NofStates[ej]);
     return true;
}


zzhang::FactorStore* zzhang::SparseEdgeNZFactor::Store(){
     zzhang::FactorStore *store = new zzhang::FactorStore(NofStates[ei] + NofStates[ej]);
     memcpy(store->data, mi, sizeof(Real) * NofStates[ei]);
     memcpy(store->data + NofStates[ei], mj, sizeof(Real) * NofStates[ej]);
     return store;
}

bool zzhang::SparseEdgeNZFactor::ReStore(zzhang::FactorStore *store)
{
     if(!store) return false;
     memcpy(mi, store->data, sizeof(Real) * NofStates[ei]);
     memcpy(mj, store->data + NofStates[ei], sizeof(Real) * NofStates[ej]);
     return true;
}

zzhang::FactorStore * zzhang::DenseEdgeFactor::Store(){
     zzhang::FactorStore *store = new FactorStore(NofStates[ei] * NofStates[ej]);
     memcpy(store->data, bij, sizeof(Real) * NofStates[ei] * NofStates[ej]);
     return store;
}

bool zzhang::DenseEdgeFactor::ReStore(zzhang::FactorStore *store)
{
     if(!store) return false;
     memcpy(bij, store->data, sizeof(Real) * NofStates[ei] * NofStates[ej]);
     return true;
}

zzhang::CGeneralSparseFactor::CGeneralSparseFactor(const std::vector<int>& Nodes,
                                                   const std::vector< std::vector<int> >& NNZs,
                                                   Real * NNZValues,
                                                   std::vector< Real * > Messages,
                                                   const std::vector< NodeFactor *> NodeFactors): m_Nodes(Nodes), m_NNZs(NNZs), m_NodeFactors(NodeFactors)
{
     this->m_NNZValues = new Real[NNZs.size()];
     memcpy(this->m_NNZValues, NNZValues, sizeof(Real) * NNZs.size());
     for(int i = 0; i < NNZs.size(); i++){
	  this->m_NNZMap[NNZs[i]] = NNZValues[i];
     }
     this->m_Messages = std::vector<Real *>(Nodes.size(), NULL);
     this->m_MaxMarginals = std::vector<Real *>(Nodes.size(), NULL);
     for(int i = 0; i < Nodes.size(); i++)
     {
	  int NofStates = m_NodeFactors[i]->m_NofStates;
	  m_Messages[i] = new Real[NofStates];
	  m_MaxMarginals[i] = new Real[NofStates];
	  if(Messages.size() > i && Messages[i] != NULL){
	       memcpy(m_Messages[i], Messages[i], sizeof(Real) * NofStates);
	  }
	  else{
	       memset(m_Messages[i], 0, sizeof(Real) * NofStates);
	  }
     }
}

zzhang::CGeneralSparseFactor::~CGeneralSparseFactor(){
     delete [] m_NNZValues;
     for(int i = 0; i < m_Messages.size(); i++) {
	  delete [] m_Messages[i];
	  delete [] m_MaxMarginals[i];
     }
}

Real zzhang::CGeneralSparseFactor::Primal(int *decode){
     std::vector<int> decodeVec(m_Nodes.size());
     double res = 0.0;
     for(int i = 0; i < m_Nodes.size(); i++){
	  decodeVec[i] = decode[m_Nodes[i]];
	  res -= m_Messages[i][decodeVec[i]];
     }
     if(m_NNZMap.find(decodeVec) != m_NNZMap.end()){
	  res += m_NNZMap[decodeVec];
     }
     return res;
}

Real zzhang::CGeneralSparseFactor::Dual(){
     return 0.0;
}

void zzhang::CGeneralSparseFactor::UpdateMessages(){
     //Particularly optimized for order 3 problems
     if(m_Messages.size() == 3){
	  int NSize = m_Messages.size();
	  std::vector< std::pair<int, Real> > MaxVs(NSize);
	  std::vector< std::pair<int, Real> > SecMaxVs(NSize);
	  std::vector< std::pair<int, Real> > ThirdVs(NSize);
	  
	  for(int i = 0; i < m_Messages.size(); i++)
	  {
	        int NofStates = m_NodeFactors[i]->m_NofStates;
		Real* bi = m_NodeFactors[i]->m_bi;
		Real* mi = m_Messages[i];
		MaxVs[i] = std::pair<int, Real>(0, bi[0] - mi[0]);
		SecMaxVs[i] = std::pair<int, Real>(0,
						   -std::numeric_limits<double>::infinity());
		ThirdVs[i] = std::pair<int, Real>(0,
						   -std::numeric_limits<double>::infinity());
		mi[0] = bi[0] - mi[0];
		for(int xi = 1; xi < NofStates; xi++)
		{
		     mi[xi] = bi[xi] - mi[xi];
		     if(mi[xi] > MaxVs[i].second){
			  ThirdVs[i] = SecMaxVs[i];
			  SecMaxVs[i] = MaxVs[i];
			  MaxVs[i] = std::pair<int, Real>(xi, mi[xi]);
		     }
		     else if(mi[xi] > SecMaxVs[i].second){
			  ThirdVs[i] = SecMaxVs[i];
			  SecMaxVs[i] = std::pair<int, Real>(xi, mi[xi]);
		     }
		     else if(mi[xi] > ThirdVs[i].second){
			  ThirdVs[i] = std::pair<int, Real>(xi, mi[xi]);
		     }
		}
	  }
	  std::vector< std::vector< std::pair<std::pair<int,int>, Real > > > LocalMaxAssign(NSize);
	  int LocalMaxIdx[3][2] = {{1,2}, {0,2}, {0,1}};
	  for(int li = 0; li < NSize; li++){
	       LocalMaxAssign[li] = std::vector< std::pair<std::pair<int,int>, Real > >(NSize);
	       std::pair<int, Real> NMaxVi = MaxVs[LocalMaxIdx[li][0]];
	       std::pair<int, Real> NSecVi = SecMaxVs[LocalMaxIdx[li][0]];

		     
	       std::pair<int, Real> NMaxVj = MaxVs[LocalMaxIdx[li][1]];
	       std::pair<int, Real> NSecVj = SecMaxVs[LocalMaxIdx[li][1]];

	       if(NMaxVi.first != NMaxVj.first){
		    LocalMaxAssign[li][0] = std::pair< std::pair<int,int>, Real>(
			 std::pair<int,int>(NMaxVi.first, NMaxVj.first),
			 NMaxVi.second + NMaxVj.second);
	       }
	       else{
		    if(NMaxVi.second - NSecVi.second < NMaxVj.second - NSecVj.second){
			 LocalMaxAssign[li][0] = std::pair< std::pair<int,int>, Real>(
			      std::pair<int,int>(NSecVi.first, NMaxVj.first),
			      NSecVi.second + NMaxVj.second);
		    }
		    else{
			 LocalMaxAssign[li][0] = std::pair< std::pair<int,int>, Real>(
			      std::pair<int,int>(NMaxVi.first, NSecVj.first),
			      NMaxVi.second + NSecVj.second);
		    }
	       }
		     
	  }
	  for(int li = 0; li < NSize; li++){
	       std::pair<int, Real> NMaxVi ;
	       std::pair<int, Real> NSecVi ;

		     
	       std::pair<int, Real> NMaxVj ;
	       std::pair<int, Real> NSecVj ;

	       int forbiddenL = LocalMaxAssign[li][0].first.first;

	       if(MaxVs[LocalMaxIdx[li][0]].first == forbiddenL){
		    NMaxVi = SecMaxVs[LocalMaxIdx[li][0]];
		    NSecVi = ThirdVs[LocalMaxIdx[li][0]];
	       }
	       else if(SecMaxVs[LocalMaxIdx[li][0]].first == forbiddenL)
	       {
		    NMaxVi = MaxVs[LocalMaxIdx[li][0]];
		    NSecVi = ThirdVs[LocalMaxIdx[li][0]];
	       }
	       else{
		    NMaxVi = MaxVs[LocalMaxIdx[li][0]];
		    NSecVi = SecMaxVs[LocalMaxIdx[li][0]];
	       }

		     
	       if(MaxVs[LocalMaxIdx[li][1]].first == forbiddenL){
		    NMaxVj = SecMaxVs[LocalMaxIdx[li][1]];
		    NSecVj = ThirdVs[LocalMaxIdx[li][1]];
	       }
	       else if(SecMaxVs[LocalMaxIdx[li][1]].first == forbiddenL)
	       {
		    NMaxVj = MaxVs[LocalMaxIdx[li][1]];
		    NSecVj = ThirdVs[LocalMaxIdx[li][1]];
	       }
	       else{
		    NMaxVj = MaxVs[LocalMaxIdx[li][1]];
		    NSecVj = SecMaxVs[LocalMaxIdx[li][1]];
	       }

		     
	       if(NMaxVi.first != NMaxVj.first){
		    LocalMaxAssign[li][1] = std::pair< std::pair<int,int>, Real>(
			 std::pair<int,int>(NMaxVi.first, NMaxVj.first),
			 NMaxVi.second + NMaxVj.second);
	       }
	       else{
		    if(NMaxVi.second - NSecVi.second < NMaxVj.second - NSecVj.second){
			 LocalMaxAssign[li][1] = std::pair< std::pair<int,int>, Real>(
			      std::pair<int,int>(NSecVi.first, NMaxVj.first),
			      NSecVi.second + NMaxVj.second);
		    }
		    else{
			 LocalMaxAssign[li][1] = std::pair< std::pair<int,int>, Real>(
			      std::pair<int,int>(NMaxVi.first, NSecVj.first),
			      NMaxVi.second + NSecVj.second);
		    }
	       }
		     
	  }
		
	  for(int li = 0; li < NSize; li++){
	       std::pair<int, Real> NMaxVi ;
	       std::pair<int, Real> NSecVi ;

		     
	       std::pair<int, Real> NMaxVj ;
	       std::pair<int, Real> NSecVj ;

	       int forbiddenL = LocalMaxAssign[li][0].first.second;

	       if(MaxVs[LocalMaxIdx[li][0]].first == forbiddenL){
		    NMaxVi = SecMaxVs[LocalMaxIdx[li][0]];
		    NSecVi = ThirdVs[LocalMaxIdx[li][0]];
	       }
	       else if(SecMaxVs[LocalMaxIdx[li][0]].first == forbiddenL)
	       {
		    NMaxVi = MaxVs[LocalMaxIdx[li][0]];
		    NSecVi = ThirdVs[LocalMaxIdx[li][0]];
	       }
	       else{
		    NMaxVi = MaxVs[LocalMaxIdx[li][0]];
		    NSecVi = SecMaxVs[LocalMaxIdx[li][0]];
	       }

		     
	       if(MaxVs[LocalMaxIdx[li][1]].first == forbiddenL){
		    NMaxVj = SecMaxVs[LocalMaxIdx[li][1]];
		    NSecVj = ThirdVs[LocalMaxIdx[li][1]];
	       }
	       else if(SecMaxVs[LocalMaxIdx[li][1]].first == forbiddenL)
	       {
		    NMaxVj = MaxVs[LocalMaxIdx[li][1]];
		    NSecVj = ThirdVs[LocalMaxIdx[li][1]];
	       }
	       else{
		    NMaxVj = MaxVs[LocalMaxIdx[li][1]];
		    NSecVj = SecMaxVs[LocalMaxIdx[li][1]];
	       }

		     
	       if(NMaxVi.first != NMaxVj.first){
		    LocalMaxAssign[li][2] = std::pair< std::pair<int,int>, Real>(
			 std::pair<int,int>(NMaxVi.first, NMaxVj.first),
			 NMaxVi.second + NMaxVj.second);
	       }
	       else{
		    if(NMaxVi.second - NSecVi.second < NMaxVj.second - NSecVj.second){
			 LocalMaxAssign[li][2] = std::pair< std::pair<int,int>, Real>(
			      std::pair<int,int>(NSecVi.first, NMaxVj.first),
			      NSecVi.second + NMaxVj.second);
		    }
		    else{
			 LocalMaxAssign[li][2] = std::pair< std::pair<int,int>, Real>(
			      std::pair<int,int>(NMaxVi.first, NSecVj.first),
			      NMaxVi.second + NSecVj.second);
		    }
	       }
		     
	  }
		
				
	  for(int i = 0; i < m_Messages.size(); i++)
	  {
	       Real *MaxMarginls = m_MaxMarginals[i];
	       Real* mi = m_Messages[i];
	       int NofStates = m_NodeFactors[i]->m_NofStates;
	       for(int xi = 0; xi < NofStates; xi++ )
	       {
		    if(xi != LocalMaxAssign[i][0].first.first &&
		       xi != LocalMaxAssign[i][0].first.second)
		    {
			 MaxMarginls[xi] = mi[xi] + LocalMaxAssign[i][0].second;
		    }
		    else if(xi == LocalMaxAssign[i][0].first.first){
			 MaxMarginls[xi] = mi[xi] + LocalMaxAssign[i][1].second;
		    }
		    else{
			 MaxMarginls[xi] = mi[xi] + LocalMaxAssign[i][2].second;
		    }
	       }
	  }
		
	  for(int i = 0; i< m_NNZs.size(); i++)
	  {
	       Real V = m_NNZValues[i];
	       for(int ni = 0; ni < m_NNZs[i].size(); ni++)
	       {
		    int xi = m_NNZs[i][ni];
		    V += m_Messages[ni][xi];
	       }
	       for(int ni = 0; ni < m_NNZs[i].size(); ni++)
	       {
		    int xi = m_NNZs[i][ni];
		    if(V > m_MaxMarginals[ni][xi]){
			 m_MaxMarginals[ni][xi] = V;
			 m_NodeFactors[ni]->m_LocalMax = xi;
		    }
	       }
	  }
	  for(int i = 0; i < m_Messages.size(); i++)
	  {
	       Real *MaxMarginls = m_MaxMarginals[i];
	       Real* mi = m_Messages[i];
	       Real* bi = m_NodeFactors[i]->m_bi;
	       int NofStates = m_NodeFactors[i]->m_NofStates;
	       for(int xi = 0; xi < NofStates; xi++ )
	       {
		    MaxMarginls[xi] /= m_Nodes.size();
		    mi[xi] = MaxMarginls[xi] - mi[xi];
		    bi[xi] = MaxMarginls[xi];
	       }
	  }
	  return;
     }


     
     std::vector<Real> MaxVs(m_Nodes.size());
     Real SumMax = 0.0;
     for(int i = 0; i < m_Messages.size(); i++){
	  int NofStates = m_NodeFactors[i]->m_NofStates;
	  Real* bi = m_NodeFactors[i]->m_bi;
	  Real* mi = m_Messages[i];
	  MaxVs[i] = -1e20;
	  for(int xi = 0; xi < NofStates; xi++){
	       mi[xi] = bi[xi] - mi[xi];
	       if(mi[xi] > MaxVs[i]){
		    MaxVs[i] = mi[xi];
		    m_NodeFactors[i]->m_LocalMax = xi;
	       }
	  }
	  SumMax += MaxVs[i];
     }
     for(int i = 0; i < m_Messages.size(); i++)
     {
	  Real base = SumMax - MaxVs[i];
	  Real *MaxMarginls = m_MaxMarginals[i];
	  Real* mi = m_Messages[i];
	  int NofStates = m_NodeFactors[i]->m_NofStates;
	  for(int xi = 0; xi < NofStates; xi++ )
	  {
	       MaxMarginls[xi] = base + mi[xi];
	  }
     }
     for(int i = 0; i< m_NNZs.size(); i++)
     {
	  Real V = m_NNZValues[i];
	  for(int ni = 0; ni < m_NNZs[i].size(); ni++)
	  {
	       int xi = m_NNZs[i][ni];
	       V += m_Messages[ni][xi];
	  }
	  for(int ni = 0; ni < m_NNZs[i].size(); ni++)
	  {
	       int xi = m_NNZs[i][ni];
	       if(V > m_MaxMarginals[ni][xi]){
		    m_MaxMarginals[ni][xi] = V;
		    m_NodeFactors[ni]->m_LocalMax = xi;
	       }
	  }
     }
     for(int i = 0; i < m_Messages.size(); i++)
     {
	  Real *MaxMarginls = m_MaxMarginals[i];
	  Real* mi = m_Messages[i];
	  Real* bi = m_NodeFactors[i]->m_bi;
	  int NofStates = m_NodeFactors[i]->m_NofStates;
	  for(int xi = 0; xi < NofStates; xi++ )
	  {
	       MaxMarginls[xi] /= m_Nodes.size();
	       mi[xi] = MaxMarginls[xi] - mi[xi];
	       bi[xi] = MaxMarginls[xi];
	  }
     }
}

void zzhang::CGeneralSparseFactor::Print(){
    
}

bool zzhang::CGeneralSparseFactor::IsGeneralFactor(){
     return false;
}
bool zzhang::CGeneralSparseFactor::GetIncludedNodes(std::vector<int> &nodes){
     nodes = m_Nodes;
     return true;
}

zzhang::FactorStore* zzhang::CGeneralSparseFactor::Store(){
     int AllSize = 0;
     for(int i = 0; i < m_NodeFactors.size(); i++)
     {
	  AllSize += m_NodeFactors[i]->m_NofStates;
     }
     zzhang::FactorStore* res = new zzhang::FactorStore(AllSize);
     int offset = 0;
     for(int i = 0; i < m_NodeFactors.size(); i++)
     {
	  memcpy(res->data + offset, m_Messages[i], sizeof(Real) * m_NodeFactors[i]->m_NofStates);
	  offset += m_NodeFactors[i]->m_NofStates;
     }
     return res;
}

bool zzhang::CGeneralSparseFactor::ReStore(zzhang::FactorStore *data){
     int offset = 0;
     for(int i = 0; i < m_NodeFactors.size(); i++)
     {
	  memcpy(m_Messages[i], data->data + offset,  sizeof(Real) * m_NodeFactors[i]->m_NofStates);
	  offset += m_NodeFactors[i]->m_NofStates;
     }
     return true;
}
