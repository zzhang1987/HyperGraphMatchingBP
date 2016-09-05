/** SubTourFactor.h --- 
 *
 * Copyright (C) 2016 Zhen Zhang
 *
 * Author: Zhen Zhang <zzhang@Zhens-MacBook-Pro.local>
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


#ifndef SUBTOURFACTOR_H
#define SUBTOURFACTOR_H 1

#include "Factor.h"

namespace zzhang{
     /**
      * A factor to forbidden subtours in TSP problem
      */
     class SubTourFactor : public CFactorBase
     {
     private:
	  /**
	   * Number of nodes in the factor
	   */
	  int NofNodes;

	  /**
	   * Number of 
	   */
	  int *NofStates;
	  /**
	   * Nodes in the facor
	   */
	  int *Nodes;

	  int *ForbiddenAssignMents;
	  /**
	   * Node beliefs, or reparametrizations
	   */
	  Real** nBeliefs;
	  /**
	   * Messages from factor to nodes
	   */
	  Real** Messages;

	  friend class CFactorGraph;

	  std::vector<NodeFactor>& m_NodeFactors;

	  
     private:
	  /**
	   * Creator
	   * @param N number of nodes
	   * @param Nodes nodes include in the factor (1xN)
	   * @param AssignMents Forbidden AssignMents
	   * @param nBeliefs, NodeBeliefs
	   */
	  SubTourFactor(int N, int *Nodes,
			int *NofStates,
			int *AssignMents,
			Real** nBeliefs,
			std::vector<NodeFactor>& NodeFactor);

	  virtual Real Primal(int* decode);
	  /**
	   * Return the dual objective. DO NOT CALL before dual updating.
	   */
	  virtual Real Dual(){
	       return 0.0;
	  }

	  virtual void UpdateMessages();

	  virtual void Print(){
	  }

	  virtual bool IsGeneralFactor(){return false;}

	  virtual bool GetIncludedNodes(std::vector<int>& nodes){
	       nodes = std::vector<int>(Nodes, Nodes+NofNodes);
	  }
	  int size(){return NofNodes;}

	  virtual ~SubTourFactor();

     private:
	  virtual zzhang::FactorStore* Store();
	  virtual bool ReStore(zzhang::FactorStore *data);
     };
}
#endif // SUBTOURFACTOR_H
