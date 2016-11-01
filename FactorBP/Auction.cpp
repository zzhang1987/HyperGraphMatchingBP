#include "Auction.h"
#include "FactorGraph.h"
void zzhang::CAuctionFactor::Auction(){
     // IsOccupied = std::vector<bool> (NofNodes, false);
	      
     for(int j = 0; j < NofNodes; j++)
     {
	  for(int i = 0; i < NofNodes; i++)
	  {
	       bi[i][j] += prices[j];
	  }
	  prices[j] = 0.0;
     }
#if 0
     for(int i = 0; i < NofNodes; i++)
     {
	  if(Evid[i] > 0)
	  {
	       int j = Evid[i] - 1;
	       if(IsOccupied[j]) return false;
	       IsOccupied[j] = true;
	  }
     }
#endif	       
     double epsilon = 5.0;
     while(epsilon > 1e-5)
     {
	  memset(AssignMent, -1, sizeof(int) * NofNodes);
	  while(1){
	       bool Finished = true;
	       for(int i = 0; i < NofNodes; i++)
	       {
		    if(AssignMent[i] == -1){
			 Finished = false;
			 break;
		    }
	       }
	       if(Finished) break;
	       auctionRound(epsilon);
	  }
	  /*
	  for(int i = 0; i < NofNodes; i++)
	  {
	       std::cout << AssignMent[i] << " ";
	  }
	  std::cout << std::endl;*/
	  double CV = m_G->ComputeObj(AssignMent);
	  for(int i = 0; i < NofNodes; i++)
	  {
	       CV -= prices[AssignMent[i]];
	  }
	  if(CV > m_G->BestDecodeV)
	  {
	       m_G->BestDecodeV = CV;
	       memcpy(m_G->m_BestDecode, AssignMent, sizeof(int) * NofNodes);
	  }
		    
	  epsilon *= 0.25;
     }

     SumPrice = 0.0;
     //  std::cout << "Prices :" << std::endl;
     for(int i = 0; i < NofNodes; i++)
     {
	  NodeFactors[i].m_LocalMax = AssignMent[i];
	  SumPrice += prices[i];
	  //std::cout << prices[i] << " ";

	  for(int j = 0; j < NofNodes; j++)
	  {
	       bi[i][j] -= prices[j];
	  }
     }
     //std::cout << std::endl;
     //return true;

}
