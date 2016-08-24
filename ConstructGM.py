import scipy.io as sio
# import cv2
import numpy as np
import numpy.matlib
from FactorGraph import *


def ConstructG(Edges1, Edges2, KP, KQ, isSparse):
    NofNodes = len(KP[0])
    NofEdges = len(Edges1[0])

    NofStates = intArray(NofNodes)
    for i in range(NofNodes):
        NofStates[i] = NofNodes

    G = CFactorGraph(NofNodes, NofStates)
    del NofStates

    G1EIdx = np.zeros([NofNodes, NofNodes])
    for i in range(NofEdges):
        ei = int(Edges1[0][i])
        ej = int(Edges1[1][i])
        G1EIdx[ei][ej] = i

    EdgeSim = KQ
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
        if (ei > ej):
            continue;
        revIdx = int(G1EIdx[ej][ei])
        EPotentials = doubleArray(NofNodes * NofNodes)
        for xij in range(NofNodes * NofNodes):
            EPotentials[xij] = 0;

        for j in range(len(Edges2[0])):
            xi = int(Edges2[0][j])
            xj = int(Edges2[1][j])
            xij = int(xi * NofNodes + xj)
            xji = int(xj * NofNodes + xi)
            EPotentials[xij] += EdgeSim[i][j]
            EPotentials[xji] += EdgeSim[revIdx][j]

            EPotentials[xi * NofNodes + xi] = -100.0
        if(isSparse):
            G.AddSparseEdgeNZ(ei, ej, EPotentials, mi, mj, NNZs, nnzIdx)
        else:
            G.AddEdge(ei, ej, EPotentials)
        del EPotentials
    if KP is not None:
        bi = doubleArray(NofNodes)
        for ni in range(NofNodes):
            for xi in range(NofNodes):
                bi[xi] = KP[ni][xi]
            G.AddNodeBelief(ni, bi)
        del bi
    G.AddAuctionFactor()
    del mi
    del mj
    return G
