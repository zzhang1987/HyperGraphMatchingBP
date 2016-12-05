import scipy.io as sio
# import cv2
import numpy as np
import numpy.matlib
from scipy.sparse import lil_matrix
from FactorBP import *
from FactorBP.FactorGraph import *

def ComputeDistAng(Edges, Pt):
    eps = 2.2204e-16
    NofE = len(Edges[0])
    EStart = Pt[Edges[0]]
    EEnd = Pt[Edges[1]]
    PtD = (EEnd - EStart).transpose()
    Dist = np.sqrt(np.sum(np.multiply(PtD, PtD), axis=0))
    angs = np.arctan(np.divide(PtD[0], PtD[1] + eps))
    return Dist, angs


def ComputeSim(Edges1, Edges2, Pt1, Pt2):
    NofE1 = len(Edges1[0])
    NofE2 = len(Edges2[0])
    simE = np.zeros([NofE1, NofE2])
    eps = 2.2204e-16

    Dist1, Angs1 = ComputeDistAng(Edges1, Pt1)
    Dist2, Angs2 = ComputeDistAng(Edges2, Pt2)

    Dst1 = np.matlib.repmat(Dist1, NofE2, 1)
    Dst2 = np.matlib.repmat(Dist2, NofE1, 1).transpose()

    Dst = np.abs(Dst1 - Dst2) / (np.max([Dst1.max(), Dst2.max()]) + eps)

    Ang1 = np.matlib.repmat(Angs1, NofE2, 1)
    Ang2 = np.matlib.repmat(Angs2, NofE1, 1).transpose()

    Ang = np.abs(Ang1 - Ang2)

    KQ = np.exp(-(Dst + Ang) / 2).transpose()
    return KQ


def ConstructDenseG(G1, Pt1, G2, Pt2):
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

            EPotentials[xi * NofNodes + xi] = -10000

        G.AddEdge(ei, ej, EPotentials)

    G.AddAuctionFactor()
    return G


def ConstructSparseG(G1, Pt1, G2, Pt2):
    Edges1 = np.nonzero(G1)
    Edges2 = np.nonzero(G2)
    NofNodes = len(G1[0])
    K = lil_matrix((NofNodes * NofNodes, NofNodes * NofNodes))
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
        if (ei > ej):
            continue;
        eij = ei * NofNodes + ej
        eji = ej * NofNodes + ei
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
            K[ei * NofNodes + xi, ej * NofNodes + xj] += EdgeSim[i][j]
            K[ei * NofNodes + xj, ej * NofNodes + xi] += EdgeSim[revIdx][j]
            EPotentials[xi * NofNodes + xi] = -10000

        G.AddSparseEdgeNZ(ei, ej, EPotentials, mi, mj, NNZs, nnzIdx)

    G.AddAuctionFactor()
    return G, K
