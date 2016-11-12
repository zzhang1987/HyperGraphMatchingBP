/** MDiverseSolver.h --- 
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


#ifndef MDIVERSESOLVER_H
#define MDIVERSESOLVER_H 1

#include "FactorGraph.h"

namespace zzhang{
     class CMdiverseSolver{
     private:
	  std::vector<CFactorGraph *> m_gphs;
	  int m_NofSols;
	  int m_NofNodes;
	  double m_Lambda;
	  std::vector< std::vector< double * > > m_AddPotentials;
     public:
	  CMdiverseSolver(int NofNodes, int NofSols, double Lambda);
	  virtual ~CMdiverseSolver();
	  void AddNodeBelief(int Nid, double *bi);
	  void UpdateMessages();
	  bool AddSparseEdgeNZ(int ei, int ej, double *data, double *mi, double *mj, int nnz, int *nnzIdx);
	  bool AddGenericGenericSparseFactor(const std::vector<int>& Nodes,
					     const std::vector< std::vector<int> >& NNZs,
					     double * NNZv);
	  void AddAuctionFactor();
	  double DualValue();
	  std::vector< std::vector<int> > GetDecode(){
	       std::vector< std::vector<int> > res(m_NofSols);
	       for(int i = 0; i < m_NofSols; i++)
	       {
		    res[i] = m_gphs[i]->GetDecode();
	       }
	       return res;
	  }
     };
}
#endif // MDIVERSESOLVER_H
