#!/usr/bin/python# FindNearSol.py --- 
# Copyright (C) 2016 Zhen Zhang
# Author: Zhen Zhang <zzhang@Zhens-MacBook-Pro.local>
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.

#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see http://www.gnu.org/licenses/.

import numpy as np
import FactorGraph as FG
import BaBSolver as BSolver
from scipy.spatial.distance import hamming
import cPickle as pickle
import MatchingGraph as MG
import FactorGraph as FG


def AdditionOfPhi(NofNodes, gamma, Phi, G):
    bi = FG.doubleArray(NofNodes)
    for i in range(NofNodes):
        for xi in range(NofNodes):
            bi[xi] = float(-Phi[i][xi] * gamma)
        G.AddNodeBelief(i, bi)

def ComputeConstraintValue(NofNodes, decode, Phi):
    res = 0;
    for i in range(NofNodes):
        res -= Phi[i][decode[i]]
    return res

def SolveConstrainedMatchingCD(NofNodes, G, Phi, X0, bestv, MaxIter=200):
    G.ResetMax()
    gamma = 0
    LastDual = 1e20;
    for iter in range(MaxIter):
        G.UpdateMessages()
        X, gamma,bestv = SolveForGamma(NofNodes, G, gamma, Phi, bestv)
        Dual = G.DualValue()
        #print('iter = %d Dual = %f Primal = %f' % (iter, Dual, bestv))
        if(np.abs(Dual - LastDual) < 1e-8):
            break;
        LastDual = Dual
        if(X is not None):
            if hamming(X,X0) > 0 :
                return X
    return None

def SolveForGamma(NofNodes, G, gamma, Phi, bestv):
    L = 0
    R = gamma + 0.5
    gamma = R
    AdditionOfPhi(NofNodes, 0.5, Phi, G)
    X = None
    while(True):
        G.RunAuction()
        cpv = G.CurrentPrimal()
        fpv = ComputeConstraintValue(NofNodes, G.GetCurrentDecode(), Phi)
        cv = cpv - R * fpv
        if(fpv > 0):
            if(cv > bestv):
                X = G.GetCurrentDecode()
                bestv = cv
            break;
        if(fpv < 0):
            AdditionOfPhi(NofNodes, R, Phi, G)
            L = R
            R = R * 2
            gamma = R
    oldgamma = gamma
    while(np.abs(L - R) > 1e-8):
        gamma = (L + R) / 2
        deltaGamma = (gamma - oldgamma)
        oldgamma = gamma
        AdditionOfPhi(NofNodes, deltaGamma, Phi, G)
        G.RunAuction()
        CX = G.GetCurrentDecode()
        CV = G.CurrentPrimal()
        fpv = ComputeConstraintValue(NofNodes, CX, Phi)
        if(fpv < 0):
            L = gamma
        else:
            R = gamma
            CV = CV - gamma * fpv
            if(CV  > bestv):
                bestv = CV
                X = CX
    return X, gamma, bestv
            

def SolveConstrainedMatching(NofNodes, G, gamma0, Phi, bestv, X0,  eps=1e-6):
    L = 0
    R = gamma0
    #bestv = G.ComputeObj(X0.tolist())
    AdditionOfPhi(NofNodes,  gamma0, Phi, G)
    X = None
    G.Solve(1000)
    while(True):
        #for i in range(MaxIter):
        #G.UpdateMessages();
        #print("")
        DualStore = G.StoreDual()
        G.ResetMax()
        res = BSolver.BaBSolver(G,100,5,0.0005,False,-1e20)
        G.ReStoreDual(DualStore)
        cpv = res.Value
        fpv = ComputeConstraintValue(NofNodes, res.Decode, Phi)

        if(np.abs(fpv) < eps):
            fpv = 0;
        if(fpv < 0):
            AdditionOfPhi(NofNodes, R, Phi, G)
            L = R
            R = R * 2
        else:
            cv = cpv - R * fpv
            CX = G.GetDecode()
            diff = np.argwhere(CX!=X0)
            NofDiff = np.prod(diff.shape)
            if(cv > bestv and NofDiff > 0):
                bestv = cv
                X = CX
                return X;
            break;
    gamma = R
    while(np.abs(L-R) > 1e-8):
        ngamma = (L + R)*1.0 / 2;
        deltagamma = (ngamma - gamma);
        gamma = ngamma
        AdditionOfPhi(NofNodes, deltagamma, Phi, G)
        G.ResetMax()
        DualStore = G.StoreDual()
        G.ResetMax()
        res = BSolver.BaBSolver(G, 100, 5, 0.000005, True, -1e20)
        G.ReStoreDual(DualStore)
        cpv = res.Value
        fpv = ComputeConstraintValue(NofNodes, res.Decode, Phi, )

        if (np.abs(fpv) < eps):
            fpv = 0;
        if(fpv < 0):
            L = gamma
        else:
            R = gamma
            CX = res.Decode;
            cv = cpv - gamma * fpv
            diff = np.argwhere(CX != X0)
            NofDiff = np.prod(diff.shape)
            if (cv > bestv and NofDiff > 0):
                bestv = cv
                X = CX
                return X;
    return X


def FindBestGamma2(NofNodes, G, X0, delta, gamma):
    Gaps = np.zeros(NofNodes);
    Xtrans = np.zeros(X0.shape, dtype=np.int32)
    for i in range(NofNodes):
        Xtrans[X0[i]] = i
    for i in range(NofNodes):
        MaxV = -1e100
        for j in range(0, NofNodes):
            if (G(j, i) > MaxV):
                MaxV = G(j, i)

        Gaps[i] = - MaxV + G(int(Xtrans[i]), i)
    nGap = np.sort(Gaps)
    deltaGamma = - nGap[delta - 1]
    for j in range(NofNodes):
        #for j in range(NofNodes):
        #    if(i == X0[j]):
        #        continue
        #    G.AddValue(j, i, -deltaGamma)
        G.AddValue(int(Xtrans[j]), j, deltaGamma)
    gamma = gamma + deltaGamma
    return gamma
        
def FindBestGamma(NofNodes, G, X0, delta, gamma):
    Gaps = np.zeros(NofNodes);
    for i in range(NofNodes):
        MaxV = -1e100
        for j in range(0, NofNodes):
            if(G(i,j) > MaxV):
                MaxV = G(i,j)
        Gaps[i] = - MaxV  + G(i,int(X0[i]))
    nGap = np.sort(Gaps)
    deltaGamma = - nGap[delta - 1]
    for i in range(NofNodes):
        #for j in range(NofNodes):
        #if(j == X0[i]):
        #    continue
        G.AddValue(i, int(X0[i]), deltaGamma)
        #G(i, X0[i]) = G(i, X0[i]) + deltaGamma
    gamma = gamma + deltaGamma
    return gamma
def FindNearMatching(NofNodes, G, delta, X0, bestv, MaxIter=1000):
    G.ResetMax()
    gamma1 = 0
    gamma2 = 0
    LastDual = 1e20;
    for iter in range(MaxIter):

        G.UpdateMessages()
        gamma1 = FindBestGamma(NofNodes, G, X0, delta, gamma1)
        gamma2 = FindBestGamma2(NofNodes, G, X0, delta, gamma2)
        G.RunAuction()

        Dual = G.DualValue() - (gamma1 + gamma2) * delta

        X = G.GetCurrentDecode()
        CV = G.CurrentPrimal()
        cpv = CV - (1 - hamming(X, X0)) * NofNodes * gamma1 + (1 - hamming(X,X0)) * NofNodes * gamma2
        if(cpv > bestv and hamming(X,X0) > 0
           and np.floor(hamming(X,X0) * NofNodes) <= delta):
            return X
        print('iter = %d Dual = %f Dis = %f, ' % (iter, Dual, hamming(X,X0) * NofNodes))
        if(np.abs(Dual - LastDual) < 1e-8):
            break;
        LastDual = Dual
    return None


def FindNearSol(NofNodes, G, X, delta, MaxIter=1000):
    #Phi = np.ones([NofNodes, NofNodes])
    G.ResetMax()
    #gamma0 = 4

    if(delta==1):
        delta = 2;


    #for i in range(NofNodes):
    #    Phi[i][X[i]] = 0

    #Phi -= delta*1.0 / NofNodes
    Xarray = FG.intArray(NofNodes)
    for i in range(NofNodes):
        Xarray[i] = int(X[i])

    cbestv = G.ComputeObj(Xarray)
    return FindNearMatching(NofNodes, G, delta, X,  cbestv)
    #return SolveConstrainedMatchingCD(NofNodes, G, Phi, X, cbestv)

    #return SolveConstrainedMatching(NofNodes,G,gamma0,Phi, cbestv, X)


def FindModes(NofNodes, G, X0, delta, MaxIter = 1000):
    X1 = None
    for iter in range(MaxIter):
        DualStore = G.StoreDual()
        X1 = FindNearSol(NofNodes, G, X0, delta)
        G.ReStoreDual(DualStore)
        if(X1 is None):
            return X0
        if(hamming(X1,X0) == 0):
            return X1
        X0 = X1
    return X1
    
    
def FindMModes(NofNodes, G, delta, N, MaxIter = 1000):
    res = dict()
    G.Solve(1000)
    DualStore = G.StoreDual()
    G.ResetMax()
    Xarray = FG.intArray(NofNodes)

    for i in range(N):
        X = np.random.permutation(NofNodes)
        X1 = FindModes(NofNodes, G, X, delta, MaxIter)
        for i in range(NofNodes):
            Xarray[i] = int(X1[i])
        X1 = np.array(X1, dtype=np.int32)
        G.ReStoreDual(DualStore)
        v = G.ComputeObj(Xarray)
        res[X1.tostring()] = v

    return res


def RunDataMModes((Fname, data, idx, NofOus)):
    car1 = data[idx]
    delta = 4
    N = 100
    LocalFeature1 = car1['features1']
    LocalFeature2 = car1['features2']
        
    PT1 = LocalFeature1[:, 0:2]
    PT2 = LocalFeature2[:, 0:2]
    
    
    orientation1 = LocalFeature1[:, 8]
    orientation2 = LocalFeature2[:, 8]
    
    GT = car1['gTruth'][0]
    
    NofInliers = len(GT)
    CMaxNofOus = np.min([LocalFeature1.shape[0], LocalFeature2.shape[0]]) - NofInliers
    CNofOus = NofOus
    if(CNofOus > CMaxNofOus):
        CNofOus = CMaxNofOus
    NofNodes = CNofOus + NofInliers
    gTruth = np.random.permutation(NofNodes)
    PT1 = PT1[gTruth, :]

    orientation1 = orientation1[gTruth]
    MG1 = MG.MatchingGraph(PT1[0:NofNodes], orientation1[0:NofNodes])
    MG2 = MG.MatchingGraph(PT2[0:NofNodes], orientation2[0:NofNodes])
    G,MFname = MG.ConstructMatchingModel(MG1, MG2, 'pas', False, True)

    MModes = FindMModes(NofNodes, G, delta, N)

    Fname = '%s_ID%d_NOus%d.pkl' % (Fname, idx, NofOus)
    f = open(Fname, "w")
    pickle.dump(MModes, f)
    pickle.dump(gTruth, f)
    pickle.dump(NofOus, f)
    f.close()
    

    
    
