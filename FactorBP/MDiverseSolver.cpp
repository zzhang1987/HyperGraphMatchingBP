/** MDiverseSolver.cpp --- 
 *
 * Copyright (C) 2016 Zhen Zhang
 *
 * Author: Zhen Zhang <zhen@server.local>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/.
 */



#include "MDiverseSolver.h"

zzhang::CMdiverseSolver::CMdiverseSolver(int NofNodes, int NofSols, double Lambda){
     m_NofSols = NofSols;
     m_NofNodes = NofNodes;
     double m_Lambda = Lambda;
     m_AddPotentials = std::vector< std::vector<double *>  > ( (NofSols) * (NofSols - 1) / 2);
     int pcnt = 0;
     int *NofStates = new int[NofNodes];
     for(int i = 0; i < NofNodes; i++)
     {
	  NofStates[i] = NofNodes;
     }
     m_gphs = std::vector<zzhang::CFactorGraph *>(m_NofSols);
     for(int i = 0; i < NofSols; i++)
     {
	  m_gphs[i] = new zzhang::CFactorGraph(m_NofNodes, NofStates);
	  for(int j = i + 1; j < NofSols; j++)
	  {
	       m_AddPotentials[pcnt] = std::vector< double * > (m_NofNodes);
	       for(int n = 0; n < m_NofNodes; n++)
	       {
		    m_AddPotentials[pcnt][n] = new double[m_NofNodes * m_NofNodes];
		    memset(m_AddPotentials[pcnt][n], 0, sizeof(double) * m_NofNodes * m_NofNodes);
		    for(int xi = 0; xi < m_NofNodes; xi++)
		    {
			 m_AddPotentials[pcnt][n][xi * m_NofNodes + xi] = - Lambda;
		    }
	       }
	       pcnt ++;
	  }
     }
     delete []NofStates;
}

zzhang::CMdiverseSolver::~CMdiverseSolver(){
     for(int i = 0; i < m_AddPotentials.size(); i++)
     {
	  for(int j = 0; j < m_AddPotentials[i].size(); j++)
	  {
	       delete []m_AddPotentials[i][j];
	  }
     }
     for(int i = 0; i < m_gphs.size(); i++)
     {
	  delete m_gphs[i];
     }
}

void zzhang::CMdiverseSolver::AddNodeBelief(int Nid, double *bi){
     for(int i = 0; i < m_gphs.size(); i++)
     {
	  m_gphs[i]->AddNodeBelief(Nid, bi);
     }
}

bool zzhang::CMdiverseSolver::AddSparseEdgeNZ(int ei, int ej,
					      double *data,
					      double *mi, double *mj,
					      int nnz, int *nnzIdx)
{
     for(int i = 0; i < m_gphs.size(); i++)
     {
	  m_gphs[i]->AddSparseEdgeNZ(ei, ej, data, mi, mj, nnz, nnzIdx);
     }
    return true;
}

bool zzhang::CMdiverseSolver::AddGenericGenericSparseFactor(const std::vector<int>& Nodes,
							    const std::vector< std::vector<int> >& NNZs,
							    double * NNZv)
{
     for(int i = 0; i < m_gphs.size(); i++)
     {
	  m_gphs[i]->AddGenericGenericSparseFactor(Nodes, NNZs, NNZv);
     }
    return true;
}

void zzhang::CMdiverseSolver::UpdateMessages(){
     int pcnt = 0;
     for(int i = 0; i < m_NofSols; i++)
     {
	  for(int j = i + 1; j < m_NofSols; j++)
	  {
	       for(int n = 0; n < m_NofNodes; n++)
	       {
		    double *bi = m_gphs[i]->m_bi[n];
		    double *bj = m_gphs[j]->m_bi[n];
		    double *bij = m_AddPotentials[pcnt][n];

		    for(int xi = 0; xi < m_NofNodes; xi++)
		    {
			 int base = xi * m_NofNodes;
			 for(int xj = 0; xj < m_NofNodes; xj++)
			 {
			      int xij = base + xj;
			      bij[xij] += bi[xi] + bj[xj];
			 }
		    }
		    for(int xi = 0; xi < m_NofNodes; xi++)
		    {
			 bi[xi] = -1e100;
			 bj[xi] = -1e100;
		    }
		    for(int xi = 0; xi < m_NofNodes; xi++)
		    {
			 int base = xi * m_NofNodes;
			 for(int xj = 0; xj < m_NofNodes; xj++)
			 {
			      int xij = base + xj;
			      if(bij[xij] > bi[xi])
				   bi[xi] = bij[xij];
			      if(bij[xij] > bj[xj])
				   bj[xj] = bij[xij];
			 }
		    }
		    for(int xi = 0; xi < m_NofNodes; xi++)
		    {
			 bi[xi] *= 0.5;
			 bj[xi] *= 0.5;
		    }
		    for(int xi = 0; xi < m_NofNodes; xi++)
		    {
			 int base = xi * m_NofNodes;
			 for(int xj = 0; xj < m_NofNodes; xj++)
			 {
			      int xij = base + xj;
			      bij[xij] -= bi[xi] + bj[xj];
			 }
		    }
	       }
	  }
     }
     for(int ni = 0; ni < m_gphs.size(); ni++)
     {
	  m_gphs[ni]->UpdateMessages();
     }
     
}


double zzhang::CMdiverseSolver::DualValue(){
     double res = 0;
     for(int i = 0; i < m_gphs.size(); i++){
	  res += m_gphs[i]->DualValue();
     }
     return res;
}

void zzhang::CMdiverseSolver::AddAuctionFactor(){
     for(int i = 0; i < m_gphs.size(); i++){
	  m_gphs[i]->AddAuctionFactor();
     }
}
