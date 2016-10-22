/** Auction.h --- 
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


#ifndef AUCTION_H
#define AUCTION_H 1
#include "PRTypes.h"
#include "Factor.h"
#include "FactorStore.h"
#include <vector>
#include <cassert>


namespace zzhang{
     class CAuctionFactor 
     {
     private:
	  int NofNodes;
	  int *NofStates;
	  Real **bi;
	  Real *prices;
	  Real SumPrice;
	  std::vector<NodeFactor>& NodeFactors;
	  std::vector<int> & Evid;
	  int *AssignMent;
	  std::vector<bool> IsOccupied;
	  friend class CFactorGraph;
	  CFactorGraph *m_G;
     CAuctionFactor(int N, int *States, Real **pi, std::vector<NodeFactor>& NFactors, std::vector<int>& evid, CFactorGraph* g): NofNodes(N),NofStates(States), bi(pi), NodeFactors(NFactors), Evid(evid)
	       {
		    m_G = g;
		    assert(N > 0);
		    prices = new Real[N];
		    AssignMent = new int[N];
		    memset(prices, 0, sizeof(Real) * N);
		    memset(AssignMent, -1, sizeof(int) * N);
	       }
	  ~CAuctionFactor(){
	       delete[] prices;
	       delete[] AssignMent;
	  }

	  FactorStore* Store(){
	       FactorStore *store = new FactorStore(NofNodes);
	       memcpy(store->data, prices, sizeof(Real) * NofNodes);
	       return store;
	  }
	  bool ReStore(FactorStore *store){
	       if(!store) return false;
	       memcpy(prices, store->data, sizeof(Real) * NofNodes);
	       SumPrice = 0;
	       for(int i = 0; i < NofNodes; i++) SumPrice += prices[i];
	       return true;
	  }

	  void Auction();

	  void auctionRound(double epsilon){
	       std::vector<int> tmpBidded;
	       std::vector<double> tmpBids;
	       std::vector<int> unAssig;
#if 0
	       for(int i = 0; i < NofNodes; i++)
	       {
		    if(Evid[i] > 0)
		    {
			 int j = Evid[i] - 1;
			 AssignMent[i] = j;
			 prices[j] = 0;
			 //std::cout << "i " << i << " j " << j << std::endl;
		    }
		    
	       }
#endif
	       for(int i = 0; i < NofNodes; i++)
	       {
		    if(AssignMent[i] == -1 /*&& Evid[i] <= 0*/)
		    {
			 unAssig.push_back(i);
			 double optValForI = -1e20;
			 double secOptValForI = -1e20;
			 int optObjForI = 0, secOptObjForI;
			 for(int j = 0; j < NofNodes; j++)
			 {
			      //if(IsOccupied[j]) continue;
			      double curVal = bi[i][j] - prices[j];
			      if(curVal > optValForI)
			      {
				   secOptValForI = optValForI;
				   secOptObjForI = optObjForI;
				   optValForI = curVal;
				   optObjForI = j;
			      }
			      else if(curVal > secOptValForI)
			      {
				   secOptValForI = curVal;
				   secOptObjForI = j;
			      }
			 }
			 double bidForI = optValForI - secOptValForI + epsilon;
			 tmpBidded.push_back(optObjForI);
			 tmpBids.push_back(bidForI);
		    }
	       }
	       
	       for(int j = 0; j <NofNodes; j++){
		    //if(IsOccupied[j]) continue;
		    std::vector<int> indices = getIndicesWithVal(&tmpBidded, j);
		    if(indices.size() != 0)
		    {
			 /* Need the highest bid for object j */
			 double highestBidForJ = tmpBids[indices[0]];
			 int i_j = indices[0];
			 for (int i = 1; i < indices.size(); i++)
			 {
			      double curVal = tmpBids.at(indices.at(i));
			      if (curVal > highestBidForJ)
			      {
				   highestBidForJ = curVal;
				   i_j = indices.at(i);
			      }
			 }
			 for(int i = 0; i < NofNodes; i++)
			 {
			      if(AssignMent[i] == j)
				   AssignMent[i] = -1;
			 }
			 AssignMent[unAssig[i_j]] = j;
			 prices[j] += highestBidForJ;
		    }
	       }
	  }

	  std::vector<int> getIndicesWithVal(std::vector<int>* v, int val){
	       std::vector<int> out;
	       for (int i = 0; i < v->size(); i++)
	       {
		    if (v->at(i) == val)
		    {
			 out.push_back(i);
		    }
	       }
	       return out;
	  }
	  
     public:
	  
     };
}
#endif // AUCTION_H
