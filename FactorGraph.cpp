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
     for(int ni = 0; ni < m_NofNodes; ni++)
     {
	  m_bi[ni] = new Real[m_NofStates[ni]];
	  memset(m_bi[ni], 0, sizeof(Real) * m_NofStates[ni]);
     }
}


void zzhang::CFactorGraph::AddNodeBelief(int Nid, Real* bi)
{
     assert(Nid < m_NofNodes);
     Real* rbi = m_bi[Nid];
     for(int xi = 0; xi < m_NofStates[Nid]; xi++)
     {
	  (*(rbi++)) += *(bi++);
     }
}

std::vector<Real> zzhang::CFactorGraph::GetBelief(int Nid)
{
     assert(Nid < m_NofNodes);
     return std::vector<Real>(m_bi[Nid],
			      m_bi[Nid] + m_NofStates[Nid]);
}
