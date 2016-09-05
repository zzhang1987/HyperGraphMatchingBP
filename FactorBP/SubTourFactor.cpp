/** SubTourFactor.cpp --- 
 *
 * Copyright (C) 2016 Zhen Zhang
 *
 * Author: Zhen Zhang <zhen@zzhang.org>
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

#include "SubTourFactor.h"

zzhang::SubTourFactor::SubTourFactor(int N, int *Nodes, int *NofStates, int *AssignMents, Real** nBeliefs, std::vector<NodeFactor>& NodeFactors): m_NodeFactors(NodeFactors)
{
     NofNodes = N;
     this->Nodes = new int[N];
     memcpy(this->Nodes, Nodes, sizeof(int) * N);
     this->NofStates = NofStates;
     this->ForbiddenAssignMents = new int[N];
     memcpy(this->ForbiddenAssignMents, AssignMents, sizeof(int) * N);
     this->nBeliefs = nBeliefs;
     Messages = new Real*[N];
     for(int i = 0; i < N; i++)
     {
	  Messages[i] = new Real[this->NofStates[this->Nodes[i]]];
	  memset(Messages[i], 0, sizeof(Real) * N);
     }
     
}

Real zzhang::SubTourFactor::Primal(int *decode){
     Real res = 0.0;
     bool IsTheSame = true;
     for(int i = 0; i < NofNodes; i++)
     {
	  if(decode[Nodes[i]] != ForbiddenAssignMents[i]) IsTheSame = false;
	  res -= Messages[i][decode[Nodes[i]]];
     }
     if(IsTheSame) return -1e20;
     else return res;
}

void zzhang::SubTourFactor::UpdateMessages(){
    

     std::vector<int> NMax(NofNodes,-1);
     std::vector<int> NSecMax(NofNodes,-1);
     std::vector<double> Gap(NofNodes, 0.0);
     std::vector<bool> IsTheSame(NofNodes,true);
     int Diff_Cnt = 0;
     double AllMax = 0.0;
     
     for(int i = 0; i < NofNodes; i++)
     {
	  int cni = Nodes[i];
	  double MaxV = -1e20;
	  double SecMaxV = -1e20;
	  for(int j = 0; j < NofStates[cni]; j++){
	       nBeliefs[cni][j] -= Messages[i][j];
	       if(nBeliefs[cni][j] > MaxV){
		    MaxV = nBeliefs[cni][j];
		    NMax[i] = j;
		    SecMaxV = MaxV;
		    NSecMax[i] = NMax[i];
	       }
	       else if(nBeliefs[cni][j] > SecMaxV){
		    SecMaxV = nBeliefs[cni][j];
		    NSecMax[i] = j;
	       }
	       
	  }
	  if(NMax[i] != ForbiddenAssignMents[i])
	  {
	       Diff_Cnt++;
	       IsTheSame[i] = false;
	       Gap[i] = 1e20;
	  }
	  else{
	       Gap[i] = MaxV - SecMaxV;
	  }
	  AllMax += MaxV;
     }
     
     Real MaxGap = -1e20;
     Real SecMaxGap = -1e20;
     int MaxIDX = -1;
    
     switch(Diff_Cnt){
     case 0:
	  for(int i = 0; i < NofNodes; i++){
	       if(Gap[i] > MaxGap){
		    MaxGap = Gap[i];
		    MaxIDX = i;
		    SecMaxGap = MaxGap;
	       }
	       else if(Gap[i] > SecMaxGap){
		    SecMaxGap = Gap[i];
	       }
	  }	  
	  for(int i = 0; i < NofNodes; i++){
	       int cni = Nodes[i];
	       double V1 = AllMax - nBeliefs[cni][NMax[i]];
	       double V2 = AllMax;
	       if(i == MaxIDX) V2 -= SecMaxGap;
	       else V2 -= MaxGap;
	       for(int xi = 0; xi < NofStates[cni]; xi++)
	       {
		    if(xi == ForbiddenAssignMents[cni])
		    {
			 double V = V2 / NofNodes;
			 Messages[i][xi] = V  - nBeliefs[cni][xi];
			 nBeliefs[cni][xi] = V ;
		    }
		    else{
			 double V = (V1 + nBeliefs[cni][xi]) / NofNodes;
			 Messages[i][xi] = V - nBeliefs[cni][xi];
			 nBeliefs[cni][xi] = V;
		    }
	       }
	       if(i==MaxIDX){
		    m_NodeFactors[cni].m_LocalMax = NSecMax[i];
	       }
	       else{
		    m_NodeFactors[cni].m_LocalMax = NMax[i];
	       }
	  }
	  
	  break;
     case 1:
	  for(int i = 0; i < NofNodes; i++){
	       if(Gap[i] > MaxGap && !IsTheSame[i]){
		    MaxGap = Gap[i];
		    MaxIDX = i;
	       }	   
	  }
	  for(int i = 0; i < NofNodes; i++)
	  {
	       int cni = Nodes[i];
	       double V1 = AllMax - nBeliefs[cni][NMax[i]];
	       for(int xi = 0; xi < NofStates[cni]; xi++)
	       {
		    double V = 0.0;
		    if(i == MaxIDX && xi == ForbiddenAssignMents[i]){
			 V = (V1 - MaxGap + nBeliefs[cni][xi]) / NofNodes;
		    }
		    else V = V1 / NofNodes;
		    Messages[i][xi] = V - nBeliefs[cni][xi];
		    nBeliefs[cni][xi] = V;
		    
	       }
	       m_NodeFactors[cni].m_LocalMax = NMax[i];
	  }
	  break;
     default:
	  for(int i = 0; i < NofNodes; i++)
	  {
	       int cni = Nodes[i];
	       double V1 = AllMax - nBeliefs[cni][NMax[i]];
	       for(int xi = 0; xi < NofStates[cni]; xi++)
	       {
		    double V = 0.0;
		    if(i == MaxIDX && xi == ForbiddenAssignMents[i]){
			 V = (V1 - MaxGap + nBeliefs[cni][xi]) / NofNodes;
		    }
		    else V = V1 / NofNodes;
		    Messages[i][xi] = V - nBeliefs[cni][xi];
		    nBeliefs[cni][xi] = V;
	       }
	       m_NodeFactors[cni].m_LocalMax = NMax[i];
	  }
     }
}


zzhang::FactorStore* zzhang::SubTourFactor::Store(){
     int TotalDim = 0;
     for(int i = 0; i < NofNodes; i++)
     {
	  int cni = Nodes[i];
	  TotalDim += NofStates[cni];
     }
     zzhang::FactorStore *store = new zzhang::FactorStore(TotalDim);
     int cnt = 0;
     for(int i = 0; i < NofNodes; i++)
     {
	  int cni = Nodes[i];
	  memcpy(store->data+cnt, Messages[i], sizeof(Real) * NofStates[cni]);
	  cnt += NofStates[cni];
     }
     return store;
}


bool zzhang::SubTourFactor::ReStore(zzhang::FactorStore *store){
     if(!store) return false;
     int cnt = 0;
     for(int i = 0; i < NofNodes; i++)
     {
	  int cni = Nodes[i];
	  memcpy(Messages[i], store->data+cnt, sizeof(Real) * NofStates[cni]);
	  cnt += NofStates[cni];
     }
     return true;
}


zzhang::SubTourFactor::~SubTourFactor(){
     delete [] Nodes;
     for(int i = 0; i < NofNodes; i++)
     {
	  delete [] Messages[i];
     }
     delete [] Messages;
}
