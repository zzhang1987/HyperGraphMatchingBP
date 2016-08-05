import scipy.io as sio
import cv2
import numpy as np
import numpy.matlib
from FactorGraph import *

def ComputeDistAng(Edges, Pt):
    eps = 2.2204e-16
    NofE = len(Edges[0])
    EStart = Pt[Edges[0]]
    EEnd = Pt[Edges[1]]
    PtD = (EEnd - EStart).transpose()
    Dist = np.sqrt(np.sum(np.multiply(PtD,PtD), axis=0))
    angs =  np.arctan(np.divide(PtD[1], PtD[0] + eps))
    return Dist, angs

def ComputeSim(Edges1, Edges2, Pt1, Pt2):
    NofE1 = len(Edges1[0])
    NofE2 = len(Edges2[0])
    simE = np.zeros([NofE1, NofE2])
    eps = 2.2204e-16
    
    Dist1,Angs1=ComputeDistAng(Edges1, Pt1)
    Dist2,Angs2=ComputeDistAng(Edges2, Pt2)

    Dst1 = np.matlib.repmat(Dist1, NofE2, 1)
    Dst2 = np.matlib.repmat(Dist2, NofE1, 1).transpose()

    Dst = np.abs(Dst1 - Dst2)/(np.max([Dst1.max(),Dst2.max()])+eps)

    Ang1 = np.matlib.repmat(Angs1, NofE2, 1)
    Ang2 = np.matlib.repmat(Angs2, NofE1, 1).transpose()

    Ang = np.abs(Ang1 - Ang2)/(np.max([Ang1.max(),Ang2.max()])+eps)

    KQ = np.exp(-(Dst+Ang)/2).transpose()
    return KQ

def ConstructSparseG(G1, Pt1, G2, Pt2):
    Edges1 = np.nonzero(G1)
    Edges2 = np.nonzero(G2)

    NofNodes = len(G1[0])
    NofStates = intArray(NofNodes);
    for i in range(NofNodes):
        NofStates[i] = NofNodes

    G1EIdx = np.zeros([NofNodes, NofNodes]);
    for i in range(len(Edges1[0])):
        ei = int(Edges1[0][i])
        ej = int(Edges1[1][i])
        G1EIdx[ei][ej] = i
    
    G = CFactorGraph(NofNodes, NofStates)
    EdgeSim = ComputeSim(Edges1, Edges2, Pt1, Pt2)
    
    NNZs = len(Edges2[0]);
    nnzIdx = intArray(2 * NNZs);
    for j in range(len(Edges2[0])):
        xi = Edges2[0][j]
        xj = Edges2[1][j]
        nnzIdx[2 * j] = int(xi)
        nnzIdx[2 * j + 1] = int(xj)

    mi = doubleArray(NofNodes);
    mj = doubleArray(NofNodes);
    for i in range(NofNodes):
        mi[i] = 0
        mj[i] = 0
    for i in range(len(Edges1[0])):
        ei = int(Edges1[0][i])
        ej = int(Edges1[1][i])
        if(ei > ej):
            continue;
        revIdx = int(G1EIdx[ej][ei])
        EPotentials = doubleArray(NofNodes * NofNodes)
        for xij in range(NofNodes*NofNodes):
            EPotentials[xij] = 0;

        for j in range(len(Edges2[0])):
            xi = int(Edges2[0][j])
            xj = int(Edges2[1][j])
            xij = int(xi * NofNodes + xj)
            xji = int(xj * NofNodes + xi)
            EPotentials[xij] += EdgeSim[i][j]
            EPotentials[xji] += EdgeSim[revIdx][j]

            EPotentials[xi * NofNodes + xi] = -10000
            
        G.AddSparseEdgeNZ(ei, ej, EPotentials, mi, mj, NNZs, nnzIdx)

    G.AddAuctionFactor()
    return G,EdgeSim
        
def RunExp(G1, Pt1, G2, Pt2):
    Edges1 = np.nonzero(G1)
    Edges2 = np.nonzero(G2)

    NofNodes = len(G1[0])
    NofStates = intArray(NofNodes);
    for i in range(NofNodes):
        NofStates[i] = NofNodes

    G1EIdx = np.zeros([NofNodes, NofNodes]);

    for i in range(len(Edges1[0])):
        ei = int(Edges1[0][i])
        ej = int(Edges1[1][i])
        G1EIdx[ei][ej] = i
    
    
    G = CFactorGraph(NofNodes, NofStates)
    G1 = CFactorGraph(NofNodes, NofStates)

    EdgeSim = ComputeSim(Edges1, Edges2, Pt1, Pt2)

    

    NNZs = len(Edges2[0]);
    nnzIdx = intArray(2 * NNZs);
    for j in range(len(Edges2[0])):
        xi = Edges2[0][j]
        xj = Edges2[1][j]
        nnzIdx[2 * j] = int(xi)
        nnzIdx[2 * j + 1] = int(xj)

    mi = doubleArray(NofNodes);
    mj = doubleArray(NofNodes);

    for i in range(NofNodes):
        mi[i] = 0
        mj[i] = 0

    for i in range(len(Edges1[0])):
        ei = int(Edges1[0][i])
        ej = int(Edges1[1][i])
        if(ei > ej):
            continue;

        revIdx = int(G1EIdx[ej][ei])
        EPotentials = doubleArray(NofNodes * NofNodes)
        for xij in range(NofNodes*NofNodes):
            EPotentials[xij] = 0;
    
        for j in range(len(Edges2[0])):
            xi = int(Edges2[0][j])
            xj = int(Edges2[1][j])
            xij = int(xi * NofNodes + xj)
            xji = int(xj * NofNodes + xi)
            EPotentials[xij] += EdgeSim[i][j]
            EPotentials[xji] += EdgeSim[revIdx][j]

            EPotentials[xi * NofNodes + xi] = -10000
            G.AddSparseEdgeNZ(ei, ej, EPotentials, mi, mj, NNZs, nnzIdx)
            G1.AddEdge(ei,ej,EPotentials)

    G.AddAuctionFactor()
    G1.AddAuctionFactor()

    #for iter in range(100):
    #    G.UpdateMessages();

    G.Solve(100)
    G1.Solve(100)
    print(G.GetDecode())

ErrorRate = np.zeros(10)

for idx in range(1):
    cnt = 0;
    SumErrorRate = 0.0;
    idx1base = (idx ) * 10;

    for d1 in range(2):
        for d2 in range(d1+1, 2):
            
            data1 = sio.loadmat('./data_chrct/' + str(idx1base+d1+1) +'.mat');
            cnt += 1
            G1 = np.array(data1['G'])
            I1 = np.array(data1['I'])
            Pt1 = np.array(data1['Pt'])
            
            data2 = sio.loadmat('./data_chrct/' + str(idx1base+d2+1) + '.mat');
            
            G2 = np.array(data2['G'])
            I2 = np.array(data2['I'])
            Pt2 = np.array(data2['Pt'])
            
            
            G,EdgeSim=ConstructSparseG(G1, Pt1, G2, Pt2)
            #G.SetVerbost(True)
            G.Solve(100)
            decode = G.GetDecode()
            ErrAssign = 0.0;
            for xi in range(len(G1)):
                if(decode[xi] != xi):
                    ErrAssign+=1;

            ErrorRate = ErrAssign / len(G1);
            SumErrorRate += ErrorRate
            #print(G.GetDecode())
    print("Char", idx, " Error Rate ", 1 - SumErrorRate/cnt)
       


