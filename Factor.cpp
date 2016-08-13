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
