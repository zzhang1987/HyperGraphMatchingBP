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

def SolveConstrainedMatching(NofNodes, G, gamma0, Phi, bestv, X0,  eps=1e-6):
    L = 0
    R = gamma0
    #bestv = G.ComputeObj(X0.tolist())
    AdditionOfPhi(NofNodes,  gamma0, Phi, G)
    X = None
    while(True):
        #for i in range(MaxIter):
        #G.UpdateMessages();
        #print("")
        DualStore = G.StoreDual()
        G.ResetMax()
        res = BSolver.BaBSolver(G,2000,5,0.005,True,-1e20)
        G.ReStoreDual(DualStore)
        cpv = res.Value
        fpv = ComputeConstraintValue(NofNodes, G.GetDecode(), Phi)

        if(np.abs(fpv) < eps):
            fpv = 0;
        if(fpv < 0):
            R = R * 2
            L = R
            AdditionOfPhi(NofNodes, R, Phi, G)
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
    while(np.abs(L-R) > 1e-6):
        ngamma = (L + R)*1.0 / 2;
        deltagamma = (ngamma - gamma);
        gamma = ngamma
        AdditionOfPhi(NofNodes, deltagamma, Phi, G)
        G.ResetMax()
        DualStore = G.StoreDual()
        G.ResetMax()
        res = BSolver.BaBSolver(G, 2000, 5, 0.005, True, -1e20)
        G.ReStoreDual(DualStore)
        cpv = res.Value
        fpv = ComputeConstraintValue(NofNodes, G.GetDecode(), Phi, )

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

                
        

def FindNearSol(NofNodes, G, X, delta, MaxIter=1000):
    Phi = np.ones([NofNodes, NofNodes])
    G.ResetMax()
    gamma0 = 4

    if(delta==1):
        delta = 2;


    for i in range(NofNodes):
        Phi[i][X[i]] = 0

    Phi -= delta*1.0 / NofNodes
    Xarray = FG.intArray(NofNodes)
    for i in range(NofNodes):
        Xarray[i] = X[i]

    cbestv = G.ComputeObj(Xarray)

    return SolveConstrainedMatching(NofNodes,G,gamma0,Phi, cbestv, X)

    
    
    
