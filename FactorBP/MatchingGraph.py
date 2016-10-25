from scipy.spatial import Delaunay
import numpy as np
from FactorBP.FactorGraph import *
from sklearn.neighbors import KDTree
import scipy.io as sio

def PermunateTriplets(T):
    T2 = np.copy(T)
    T3 = np.copy(T)
    T4 = np.copy(T)
    T5 = np.copy(T)
    T6 = np.copy(T)
    
    T2[:,0] = T[:,1]
    T2[:,1] = T[:,0]
    
    T3[:,0] = T[:,1]
    T3[:,1] = T[:,2]
    T3[:,2] = T[:,0]
    
    T4[:,0] = T[:,2]
    T4[:,1] = T[:,1]
    T4[:,2] = T[:,0]
    
    T5[:,0] = T[:,2]
    T5[:,1] = T[:,0]
    T5[:,2] = T[:,1]
    
    T6[:,0] = T[:,0]
    T6[:,1] = T[:,2]
    T6[:,2] = T[:,1]
    
    res = np.append(T,T2,axis = 0)
    res = np.append(res,T3,axis = 0)
    res = np.append(res,T4,axis = 0)
    res = np.append(res,T5,axis = 0)
    res = np.append(res,T6,axis = 0)
    
    return res

def computeTripletsFeatureSinAlpha(Points, T):
    vecX = np.zeros(3)
    vecY = np.zeros(3)
    F = np.zeros(3)
    if ((T[0] == T[1]) or (T[0] == T[2]) or (T[1] == T[2])):
        F = -10 * np.ones(3)
        return F
    for idx in range(3):
        vecX[idx] = Points[T[(idx + 1) % 3]][0] - Points[T[idx]][0]
        vecY[idx] = Points[T[(idx + 1) % 3]][1] - Points[T[idx]][1]
        length = np.linalg.norm([vecX[idx], vecY[idx]])
        if (length != 0):
            vecX[idx] /= length
            vecY[idx] /= length
        else:
            vecX[idx] = 0
            vecY[idx] = 0

    for idx in range(3):
        F[idx] = vecX[((idx+1)%3)]*vecY[idx]-vecY[((idx+1)%3)]*vecX[idx]

    return F
def computeTripletsFeatureSimple(Points, T):
    vecX = np.zeros(3)
    vecY = np.zeros(3)
    F = np.zeros(3)
    if((T[0]==T[1]) or (T[0]==T[2]) or (T[1]==T[2])):
        F = -10 * np.ones(3)
        return F
    for idx in range(3):
        vecX[idx] = Points[T[(idx + 1)%3]][0] - Points[T[idx]][0]
        vecY[idx] = Points[T[(idx + 1)%3]][1] - Points[T[idx]][1]
        length = np.linalg.norm([vecX[idx], vecY[idx]])
        if(length != 0):
            vecX[idx] /= length
            vecY[idx] /= length
        else:
            vecX[idx] = 0
            vecY[idx] = 0
            
            
    for idx in range(3):
        F[idx] = np.arctan2(vecX[idx], vecY[idx])
        
    return F
def computeEdgeFeatureSimple(Points, E1, E2):
    Vec = Points[E1] - Points[E2]
    F = np.zeros(4)
    F[0] = np.linalg.norm(Vec)
    F[2] = F[0]
    F[1] = np.arctan2(Vec[0], Vec[1])
    F[3] = np.arctan2(-Vec[0], Vec[1])
    return F


def ComputeFeatureDistance(F1, F2):
    res = np.zeros([F1.shape[0], F2.shape[0]])
    for i in range(F1.shape[0]):
        for j in range(F2.shape[0]):
            res[i][j] = (np.linalg.norm(F1[i] - F2[j]))
    return res
def ComputeMultiAngleDistance(F1,F2):
    res = np.zeros([F1.shape[0], F2.shape[0]])
    for i in range(F1.shape[0]):
        for j in range(F2.shape[0]):
            AngD = F1[i] - F2[j]
            for ai in range(AngD.shape[0]):
                if(AngD[ai] < 0):
                    AngD[ai] = AngD[ai]
                if(AngD[ai] > np.pi):
                    AngD[ai] = 2 * np.pi - AngD[ai]
            res[i][j] = np.linalg.norm(AngD)
    return res
def ComputeAngleDistance(F1, F2):
    res = np.zeros([F1.shape[0], F2.shape[0]])
    for i in range(F1.shape[0]):
        for j in range(F2.shape[0]):
            res[i][j] = (F1[i] - F2[j])
            if(res[i][j] < 0):
                res[i][j] = -res[i][j]
            if(res[i][j] > np.pi):
                res[i][j] -= np.pi
    return res
def ComputeKQ(G1, G2, Type):
    NofEdges1 = G1.Edges.shape[0]
    NofEdges2 = G2.Edges.shape[0] * 2
    KQ = np.zeros([NofEdges1, NofEdges2])
    if(Type == 'pas'):
        distTable = ComputeFeatureDistance(G1.EdgeFeature[:, 0],
                                           np.append(G2.EdgeFeature[:,0],G2.EdgeFeature[:,2]))
        for i in range(NofEdges1):
            for j in range(G2.Edges.shape[0]):
                distTable[i][j] /= (np.min([G1.EdgeFeature[i][0], G2.EdgeFeature[j][0]]) + 1e-6)
                distTable[i][G2.Edges.shape[0] + j] /= (np.min([G1.EdgeFeature[i][0], G2.EdgeFeature[j][2]]) + 1e-6)


        agdistTable = ComputeAngleDistance(G1.EdgeFeature[:, 1], np.append(G2.EdgeFeature[:,1], G2.EdgeFeature[:,3]))


        agdistTable1 = ComputeAngleDistance(G1.EdgeFeature[:, 3], np.append(G2.EdgeFeature[:,3], G2.EdgeFeature[:,1]))

        #MinDis = np.mean([np.min(G1.EdgeFeature[:,0]), np.min(G2.EdgeFeature[:,0])])
        #distTable /= (MinDis + 1e-6)

        KQ = np.exp(-(distTable + agdistTable)/2)
        KQ1 = np.exp(-(distTable + agdistTable1)/2)
        KQ = KQ + KQ1
    if(Type == 'cmu'):
        distTable = ComputeFeatureDistance(G1.EdgeFeature[:, 0],
                                           np.append(G2.EdgeFeature[:,0],G2.EdgeFeature[:,2]))
        
        KQ = np.exp(-(distTable**2)/2500) * 2
        
    return KQ

def ComputeKT(G1,G2):
    distTable = ComputeMultiAngleDistance(G1.TFeature, G2.PTFeature)
    KT = np.exp(-distTable/np.mean(distTable))
    return KT

def ConstructMatchingModelRandom(G1, G2, Type, AddTriplet):
    KP = ComputeFeatureDistance(G1.PFeature, G2.PFeature)
    KQ = ComputeKQ(G1, G2, Type)

    NP1 = G1.NofNodes
    NP2 = G2.NofNodes
    nT = np.floor(NP1 * NP2)
    t1 = np.floor(np.random.rand(3, nT) * NP1)
    while (True):
        probFound = False;
        for i in range(3):
            ind = (t1[i, :] == t1[(i + 1) % 3, :])
            if (np.sum(ind) != 0):
                idxs = np.nonzero(ind)
                t1[i][idxs] = np.floor(np.random.rand(1, len(idxs[0])) * NP1);
                probFound = True;
        if (probFound == False):
            break;

    t1 = t1.transpose()
    T = np.sort(t1,axis = 1)
    T = T[np.lexsort(np.fliplr(T).T)]
    NRepeated = np.ones(T.shape[0], dtype=int)
    for i in range(1,T.shape[0]):
        if(np.sum(np.abs(T[i] -T[i-1])) == 0):
            NRepeated[i] = 0

    NRepTri = np.nonzero(NRepeated)
    T = T[NRepTri]
    TF1 = np.zeros([T.shape[0],3])
    for ti in range(T.shape[0]):
        TF1[ti] = computeTripletsFeatureSinAlpha(G1.P, T[ti])


    NofT2 = G2.NofNodes * (G2.NofNodes - 1) * (G2.NofNodes - 2)
    T2 = np.zeros([NofT2, 3], dtype=int)
    TF2 = np.zeros([6 * NofT2, 3], dtype=float)
    T2Cnt = 0
    for i1 in range(G2.NofNodes):
        for i2 in range(i1+1, G2.NofNodes):
            for i3 in range(i2 + 1, G2.NofNodes):
                T2[T2Cnt][0] = i1
                T2[T2Cnt][1] = i2
                T2[T2Cnt][2] = i3
                T2Cnt +=1
    T2 = PermunateTriplets(T2)

    for ti in range(T2.shape[0]):
        TF2[ti] = computeTripletsFeatureSinAlpha(G2.P, T2[ti])

    kdt = KDTree(TF2, metric='euclidean')
    nNN = T.shape[0] * 2
    [distT, indicesT] = kdt.query(TF1, k=nNN, return_distance=True)
    distT = np.exp(- (distT / np.mean(distT) ))
    KP = np.exp(-KP)
    NofNodes = G1.NofNodes
    NofStates = intArray(NofNodes)
    for i in range(NofNodes):
        NofStates[i] = NofNodes
    G = CFactorGraph(NofNodes, NofStates)
    bi = doubleArray(NofNodes)
    for ni in range(NofNodes):
        for xi in range(NofNodes):
            bi[xi] = float(KP[ni][xi])
        G.AddNodeBelief(ni, bi)
    nnzEdgeIdx = VecVecInt(KQ.shape[1])
    for ni in range(G2.Edges.shape[0]):
        CurrentAssign = VecInt(2)
        CurrentAssign[0] = int(G2.Edges[ni][0])
        CurrentAssign[1] = int(G2.Edges[ni][1])
        InvCurrentAssign = VecInt(2)
        InvCurrentAssign[0] = int(G2.Edges[ni][1])
        InvCurrentAssign[1] = int(G2.Edges[ni][0])
        nnzEdgeIdx[ni] = CurrentAssign
        nnzEdgeIdx[ni + G2.Edges.shape[0]] = InvCurrentAssign

    for ei in range(KQ.shape[0]):
        CEdgeVec = VecInt(2)
        CEdgeVec[0] = int(G1.Edges[ei][0])
        CEdgeVec[1] = int(G1.Edges[ei][1])
        CurrentNNZV = doubleArray(KQ.shape[1])
        for xij in range(KQ.shape[1]):
            CurrentNNZV[xij] = KQ[ei][xij]
        G.AddGenericGenericSparseFactor(CEdgeVec, nnzEdgeIdx, CurrentNNZV)


    for ti in range(distT.shape[0]):
        CTripletsVec = VecInt(3)
        CTripletsVec[0] = int(T[ti][0])
        CTripletsVec[1] = int(T[ti][1])
        CTripletsVec[2] = int(T[ti][2])
        nnzTripIdx = VecVecInt(distT.shape[1])
        nnzTripV = doubleArray(distT.shape[1])
        for xijk in range(distT.shape[1]):
            cIdxVec = VecInt(3)
            cIdxVec[0] = int(T2[indicesT[ti][xijk]][0]);
            cIdxVec[1] = int(T2[indicesT[ti][xijk]][1]);
            cIdxVec[2] = int(T2[indicesT[ti][xijk]][2]);
            nnzTripIdx[xijk] = cIdxVec
            nnzTripV[xijk] = 6 * distT[ti][xijk]
        G.AddGenericGenericSparseFactor(CTripletsVec, nnzTripIdx, nnzTripV)

    G.AddAuctionFactor()

    return G



def ConstructMatchingModel(G1, G2, Type, AddTriplet):
    KP = ComputeFeatureDistance(G1.PFeature, G2.PFeature)
    KQ = ComputeKQ(G1, G2, Type)
    KT = ComputeKT(G1, G2)
    KP = np.exp(-(KP))
    NofNodes = G1.NofNodes
    NofStates = intArray(NofNodes)
    for i in range(NofNodes):
        NofStates[i] = NofNodes
    G = CFactorGraph(NofNodes, NofStates)
    bi = doubleArray(NofNodes)
    for ni in range(NofNodes):
        for xi in range(NofNodes):
            bi[xi] = float(KP[ni][xi])
        G.AddNodeBelief(ni, bi)
    nnzEdgeIdx = VecVecInt(KQ.shape[1])
    for ni in range(G2.Edges.shape[0]):
        CurrentAssign = VecInt(2)
        CurrentAssign[0] = int(G2.Edges[ni][0])
        CurrentAssign[1] = int(G2.Edges[ni][1])
        InvCurrentAssign = VecInt(2)
        InvCurrentAssign[0] = int(G2.Edges[ni][1])
        InvCurrentAssign[1] = int(G2.Edges[ni][0])
        nnzEdgeIdx[ni] = CurrentAssign
        nnzEdgeIdx[ni + G2.Edges.shape[0]] = InvCurrentAssign

    nnzTripIdx = VecVecInt(KT.shape[1])
    for ni in range(KT.shape[1]):
        CurrentAssign = VecInt(3)
        CurrentAssign[0] = int(G2.PermunatedTriplets[ni][0])
        CurrentAssign[1] = int(G2.PermunatedTriplets[ni][1])
        CurrentAssign[2] = int(G2.PermunatedTriplets[ni][2])
        nnzTripIdx[ni] = CurrentAssign


    NNZs = KQ.shape[1]
    nnzIdx = intArray(2 * NNZs)
    for j in range(G2.Edges.shape[0]):
        xi = G2.Edges[j][0]
        xj = G2.Edges[j][1]
        nnzIdx[2 * j] = int(xi)
        nnzIdx[2 * j + 1] = int(xj)
    baseIdx = 2 * G2.Edges.shape[0]
    for j in range(G2.Edges.shape[0]):
        xi = G2.Edges[j][1]
        xj = G2.Edges[j][0]
        nnzIdx[baseIdx + 2 * j] = int(xi)
        nnzIdx[baseIdx + 2 * j + 1] = int(xj)
    
    mi = doubleArray(NofNodes);
    mj = doubleArray(NofNodes);
    for i in range(NofNodes):
        mi[i] = 0
        mj[i] = 0
    for ei in range(KQ.shape[0]):
        EPotentials = doubleArray(NofNodes * NofNodes)
        for xij in range(NofNodes * NofNodes):
            EPotentials[xij] = 0;
        cei = G1.Edges[ei][0]
        cej = G1.Edges[ei][1]
        CEdgeVec = VecInt(2)
        CEdgeVec[0] = int(G1.Edges[ei][0])
        CEdgeVec[1] = int(G1.Edges[ei][1])
        CurrentNNZV = doubleArray(KQ.shape[1])
        for xij in range(KQ.shape[1]):
            rxij = int(nnzIdx[2 * xij] * NofNodes + nnzIdx[2 * xij + 1])
            EPotentials[rxij] = KQ[ei][xij]
            CurrentNNZV[xij] = KQ[ei][xij]
        #G.AddGenericGenericSparseFactor(CEdgeVec, nnzEdgeIdx, CurrentNNZV)
        G.AddSparseEdgeNZ(cei, cej, EPotentials, mi, mj, NNZs, nnzIdx)


    if(AddTriplet == False):
        G.AddAuctionFactor()
        return G

    for ti in range(KT.shape[0]):
    #for ti in range(0):
        CTripletsVec = VecInt(3)
        CTripletsVec[0] = int(G1.Triplets[ti][0])
        CTripletsVec[1] = int(G1.Triplets[ti][1])
        CTripletsVec[2] = int(G1.Triplets[ti][2])
        CurrentNNZV = doubleArray(KT.shape[1])
        for xijk in range(KT.shape[1]):
            CurrentNNZV[xijk] = 0.32 * KT[ti][xijk]
        G.AddGenericGenericSparseFactor(CTripletsVec, nnzTripIdx, CurrentNNZV)

    G.AddAuctionFactor()

    MatRes = {}
    MatRes['Triplets'] = G1.Triplets
    MatRes['NTriplets'] = G2.PermunatedTriplets
    MatRes['Similarity'] = KT
    MatRes['Edges'] = G1.Edges
    MatRes['NEdges'] = G2.Edges
    MatRes['KQ'] = KQ
    MatRes['KP'] = KP
    MatRes['GT'] = range(NofNodes)
    sio.savemat('Temp.mat', MatRes)

    return G;



class MatchingGraph:
    def __init__(self, P, PFeature):
        self.P = P
        self.PFeature = PFeature
        self.NofNodes = P.shape[0]
        self.Delaunay()
        self.AddEdges()
        self.ComputeFeature()
    def ComputeFeature(self):
        self.ComputeEdgeFeature()
        self.ComputeTripletFeature()
    def ComputeEdgeFeature(self):
        self.EdgeFeature = np.zeros([self.Edges.shape[0], 4])
        for ei in range(self.Edges.shape[0]):
            self.EdgeFeature[ei] = computeEdgeFeatureSimple(self.P, self.Edges[ei][0], self.Edges[ei][1])
    def ComputeTripletFeature(self):
        self.TFeature = np.zeros([self.Triplets.shape[0], 3])
        for ti in range(self.Triplets.shape[0]):
            self.TFeature[ti] = computeTripletsFeatureSimple(self.P, self.Triplets[ti])
                                     
        self.PTFeature = np.zeros([self.PermunatedTriplets.shape[0], 3])
        for ti in range(self.PermunatedTriplets.shape[0]):
            self.PTFeature[ti] = computeTripletsFeatureSimple(self.P, self.PermunatedTriplets[ti])
    def AddEdges(self):
        G = np.zeros([self.NofNodes, self.NofNodes])
        for i in range(self.Triplets.shape[0]):
            G[self.Triplets[i][0]][self.Triplets[i][1]] = 1
            G[self.Triplets[i][1]][self.Triplets[i][0]] = 1
            G[self.Triplets[i][0]][self.Triplets[i][2]] = 1
            G[self.Triplets[i][2]][self.Triplets[i][0]] = 1
            G[self.Triplets[i][1]][self.Triplets[i][2]] = 1
            G[self.Triplets[i][2]][self.Triplets[i][1]] = 1
        NofEdges = int(np.sum(G)/2)
        self.Edges = np.zeros([NofEdges, 2], dtype=int)
        Ecnt = 0
        for i in range(self.NofNodes):
            for j in range(i+1, self.NofNodes):
                if(G[i][j] == 1):
                    self.Edges[Ecnt][0] = i
                    self.Edges[Ecnt][1] = j
                    Ecnt += 1
        
                
    def Delaunay(self):
        tri1 = Delaunay(self.P)
        self.Triplets = tri1.simplices
        self.PermunatedTriplets = PermunateTriplets(self.Triplets)
        
