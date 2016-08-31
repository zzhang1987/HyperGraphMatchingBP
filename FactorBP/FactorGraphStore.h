/** FactorGraphStore.h --- 
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


#ifndef FACTORGRAPHSTORE_H
#define FACTORGRAPHSTORE_H 1
#include <vector>
#include "FactorStore.h"
namespace zzhang{
     class FactorGraphDualStore{
     public:
	  std::vector<FactorStore *> NodeStores;
	  std::vector<FactorStore *> FactorStores;
	  FactorStore *AuctionStore;
	  std::vector<int> Evid;
	  FactorGraphDualStore(int NofNodes, int NofFactors){
	       assert(NofNodes > 0 && NofFactors >= 0);
	       NodeStores = std::vector<FactorStore *>(NofNodes, NULL);
	       FactorStores = std::vector<FactorStore *>(NofFactors, NULL);
	       AuctionStore = NULL;
	  }
	  virtual ~FactorGraphDualStore(){
	       //std::cout << "Deleted" << std::endl;
	       for(int i = 0; i < NodeStores.size(); i++)
	       {
		    if(NodeStores[i])
			 delete NodeStores[i];
	       }
	       for(int i = 0; i < FactorStores.size(); i++){
		    if(FactorStores[i])
			 delete FactorStores[i];
	       }
	       if(AuctionStore) delete AuctionStore;
	  }
     };
}
#endif // FACTORGRAPHSTORE_H
