/** testMain.cpp --- 
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



#include "FactorGraph.h"
#include <boost/shared_array.hpp>
#include <iostream>
int main(int argc, char **argv){
    int NofNodes = 3;
    int NofStates[3] = {3,3,3};
    double NodeBliefs[3][3] = {{0.05,0.43,0.23}, {0.124,0.179683,0.169873}, {0.293847,0.91687,0.103789} };
    double NNZV[2]={0.05,0.45};
    zzhang::CFactorGraph FG(NofNodes, NofStates);
    FG.AddNodeBelief(0, NodeBliefs[0]);
    FG.AddNodeBelief(1, NodeBliefs[1]);
    FG.AddNodeBelief(2, NodeBliefs[2]);
    std::vector< std::vector<int> > NNZs(2);
    std::vector<int> Nodes(3);
    Nodes[0] = 0;
    Nodes[1] = 1;
    Nodes[2] = 2;
    NNZs[0] = std::vector<int>(3);
    NNZs[0][0] = 0;
    NNZs[0][1] = 1;
    NNZs[0][2] = 2;
    NNZs[1] = std::vector<int>(3);
    NNZs[1][0] = 1;
    NNZs[1][1] = 2;
    NNZs[1][1] = 0;
    
    FG.AddGenericGenericSparseFactor(Nodes, NNZs, NNZV);
    
    FG.AddAuctionFactor();
    FG.SetVerbose(true);
    FG.SetMinDualDecrease(1e-8);
    FG.Solve(100);
    
    
    
    std::cout << "Test " << std::endl;
}
