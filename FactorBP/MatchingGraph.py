from scipy.spatial import Delaunay
import numpy as np
from FactorGraph import CFactorGraph, intArray, doubleArray, VecInt, VecVecInt
from sklearn.neighbors import KDTree
import scipy.io as sio
import tempfile as tmp
from scipy.spatial.distance import hamming
def rand_rotation_matrix(deflection=1.0, randnums=None):
    """
    Creates a random rotation matrix.
    
    deflection: the magnitude of the rotation. For 0, no rotation; for 1, competely random
    rotation. Small deflection => small perturbation.
    randnums: 3 random numbers in the range [0, 1]. If `None`, they will be auto-generated.
    """
    # from http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
    
    if randnums is None:
        randnums = np.random.uniform(size=(3,))
        
    theta, phi, z = randnums
  
    theta = theta * 2.0*deflection*np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0*np.pi  # For direction of pole deflection.
    z = z * 2.0*deflection  # For magnitude of pole deflection.
    
    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.
    
    r = np.sqrt(z)
    Vx, Vy, Vz = V = (
        np.sin(phi) * r,
        np.cos(phi) * r,
        np.sqrt(2.0 - z)
        )
    
    st = np.sin(theta)
    ct = np.cos(theta)
    
    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))
    
    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
    
    M = (np.outer(V, V) - np.eye(3)).dot(R)
    return M



def GenRandomMatchingPoints(NofInliers, Scale,  Noise, NofOutliers, theta = 0):
    MaxSize = 1000
    PT1 = np.random.rand(NofInliers, 2) * MaxSize
    #PT1[:,2] *= 0

    #PT1Homo = np.append(PT1, np.ones([NofInliers, 1]), axis = 1)

    PT1Homo = PT1.transpose()

    
    TransMat3 = np.zeros([2, 2])
    TransMat3[0][0] = np.cos(theta)
    TransMat3[0][1] = -np.sin(theta)
    TransMat3[1][0] = np.sin(theta)
    TransMat3[1][1] = np.cos(theta)


    

    TransMat = Scale * TransMat3
    #TransMat = rand_rotation_matrix()
    #TransMat[2][2] = 1

    #TransMat[0][0] = np.cos(theta) * Scale
    #TransMat[0][1] = -np.sin(theta) * Scale
    #TransMat[1][0] = np.sin(theta) * Scale
    #TransMat[1][1] = np.cos(theta) * Scale
    print(TransMat)
    PT2Trans = TransMat.dot(PT1Homo) + np.random.normal(0, Noise, [2, NofInliers]) #Noise * np.random.rand(2, NofInliers)
    PT2Homo = PT2Trans.transpose()
    #PT2 = PT2Homo[:,0:2]

    PT2 = PT2Homo


    Ous1 = np.random.rand(NofOutliers, 2) * MaxSize * 4
    Ous2 = np.random.rand(NofOutliers, 2) * MaxSize * 4 * Scale
    #PT1 = PT1[:,0:2]
    #PT1[:,0] *= 0.8
    PT11 = np.append(PT1, Ous1, axis = 0)
    PT22 = np.append(PT2, Ous2, axis = 0)
    
    return PT11,PT22

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


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def computeTripletsFeatureSimple(Points, T):
    vecX = np.zeros(3)
    vecY = np.zeros(3)
    F = np.zeros(3)
    if((T[0]==T[1]) or (T[0]==T[2]) or (T[1]==T[2])):
        F = -10 * np.ones(3)
        return F
    for idx in range(3):
        P1 = Points[T[(idx + 1)%3]] -  Points[T[idx]]
        P2 = Points[T[(idx + 2)%3]] - Points[T[idx]]
  
    
        F[idx] = angle_between(P1,P2)
        # length = np.linalg.norm([vecX[idx], vecY[idx]])
        # if(length != 0):
        #    vecX[idx] /= length
        #    vecY[idx] /= length
        # else:
        #    vecX[idx] = 0
        #    vecY[idx] = 0
    # for idx in range(3):
        # F[idx] = np.arctan2(vecX[idx], vecY[idx])
        
    return F


def computeEdgeFeatureSimple(Points, E1, E2):
    Vec = Points[E1] - Points[E2]
    F = np.zeros(4)
    F[0] = np.linalg.norm(Vec)
    F[2] = F[0]
    F[1] = np.arctan2(Vec[0], Vec[1])
    F[3] = np.arctan2(-Vec[0], -Vec[1])
    return F


def ComputeFeatureDistance(F1, F2, dis='L2'):
    res = np.zeros([F1.shape[0], F2.shape[0]])
    for i in range(F1.shape[0]):
        for j in range(F2.shape[0]):
            if(dis == 'L2'):
                res[i][j] = (np.linalg.norm(F1[i] - F2[j]))
            elif(dis == 'hamming'):
                res[i][j] = hamming(F1[i],F2[j])
    return res
def ComputeMultiAngleDistance(F1,F2):
    res = np.zeros([F1.shape[0], F2.shape[0]])
    for i in range(F1.shape[0]):
        for j in range(F2.shape[0]):
            AngD = F1[i] - F2[j]
            for ai in range(AngD.shape[0]):
                if(AngD[ai] < 0):
                    AngD[ai] = -AngD[ai]
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
                res[i][j] = 2 * np.pi - res[i][j]
    return res

def ComputeKQ(G1, G2, Type):
    NofEdges1 = G1.Edges.shape[0]
    NofEdges2 = G2.Edges.shape[0] * 2
    KQ = np.zeros([NofEdges1, NofEdges2])
    if(Type == 'pasDis'):
        distTable = ComputeFeatureDistance(G1.EdgeFeature[:, 0],
                                           np.append(G2.EdgeFeature[:,0],G2.EdgeFeature[:,2]))
        for i in range(NofEdges1):
            for j in range(G2.Edges.shape[0]):
                distTable[i][j] /= (np.min([G1.EdgeFeature[i][0], G2.EdgeFeature[j][0]]) + 1e-6)
                distTable[i][G2.Edges.shape[0] + j] /= (np.min([G1.EdgeFeature[i][0], G2.EdgeFeature[j][2]]) + 1e-6)
        agdistTable = ComputeAngleDistance(G1.EdgeFeature[:, 1], np.append(G2.EdgeFeature[:,1], G2.EdgeFeature[:,3]))

        KQ = np.exp(-(distTable + agdistTable)/2) * 2
    if(Type == 'pasDisOnly'):
        distTable = ComputeFeatureDistance(G1.EdgeFeature[:, 0],
                                           np.append(G2.EdgeFeature[:,0],G2.EdgeFeature[:,2]))
        for i in range(NofEdges1):
            for j in range(G2.Edges.shape[0]):
                distTable[i][j] /= (np.min([G1.EdgeFeature[i][0], G2.EdgeFeature[j][0]]) + 1e-6)
                distTable[i][G2.Edges.shape[0] + j] /= (np.min([G1.EdgeFeature[i][0], G2.EdgeFeature[j][2]]) + 1e-6)
        KQ = np.exp(-distTable) * 2
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
        
    if(Type == 'syn'):
        distTable = ComputeFeatureDistance(G1.EdgeFeature[:, 0],
                                           np.append(G2.EdgeFeature[:,0],G2.EdgeFeature[:,2]))
        for i in range(NofEdges1):
            for j in range(G2.Edges.shape[0]):
                distTable[i][j] /= (np.min([G1.EdgeFeature[i][0], G2.EdgeFeature[j][0]]) + 1e-6)
                distTable[i][G2.Edges.shape[0] + j] /= (np.min([G1.EdgeFeature[i][0], G2.EdgeFeature[j][2]]) + 1e-6)
        KQ = np.exp(-(distTable)) * 2
        #KQ = np.zeros(distTable.shape)
    if(Type == 'syn_noedge'):
        distTable = ComputeFeatureDistance(G1.EdgeFeature[:, 0],
                                           np.append(G2.EdgeFeature[:,0],G2.EdgeFeature[:,2]))
        for i in range(NofEdges1):
            for j in range(G2.Edges.shape[0]):
                distTable[i][j] /= (np.min([G1.EdgeFeature[i][0], G2.EdgeFeature[j][0]]) + 1e-6)
                distTable[i][G2.Edges.shape[0] + j] /= (np.min([G1.EdgeFeature[i][0], G2.EdgeFeature[j][2]]) + 1e-6)
        KQ = np.exp(-(distTable)) * 2
        KQ = np.ones(distTable.shape)
    if (Type == 'topo'): # add by Lee at 13:44PM 9th November
        distTable = ComputeFeatureDistance(G1.EdgeFeature[:, 0],
                                           np.append(G2.EdgeFeature[:, 0], G2.EdgeFeature[:, 2]))
        #for i in range(NofEdges1):
        #    for j in range(G2.Edges.shape[0]):
        #        distTable[i][j] /= (np.min([G1.EdgeFeature[i][0], G2.EdgeFeature[j][0]]) + 1e-6)
        #        distTable[i][G2.Edges.shape[0] + j] /= (np.min([G1.EdgeFeature[i][0], G2.EdgeFeature[j][2]]) + 1e-6)
        # KQ = np.exp(-(distTable)) * 2
        KQ = np.ones(distTable.shape)
    return KQ

def ComputeKT(G1,G2):
    distTable = ComputeMultiAngleDistance(G1.TFeature, G2.PTFeature)
    KT = np.exp(-distTable/np.mean(distTable) )
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
    nNN = T.shape[0]
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

def ComputeSimilarity(G1, G2, Type):
    if(Type == 'syn'):
        KP = ComputeFeatureDistance(G1.PFeature, G2.PFeature, dis='hamming')
    else:
        KP = ComputeFeatureDistance(G1.PFeature, G2.PFeature)

    KQ = ComputeKQ(G1, G2, Type)
    KT = 2 * ComputeKT(G1, G2)
    KP = np.exp(-(KP))
    return KP,KQ,KT

# Modifyed by Lee at 12:34PM 4th November
def ConstructMatchingModel(G1, G2, Type, AddTriplet = True, AddEdge = True):
    KP,KQ,KT = ComputeSimilarity(G1,G2,Type)
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
    # NNZEdge & NNZTrip?
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
    # Save Matching Model as Temp.mat
    MatRes = {}
    TmpFName = tmp.mktemp(suffix='.mat')
    if AddTriplet:
        MatRes['Triplets'] = G1.Triplets
        MatRes['NTriplets'] = G2.PermunatedTriplets
        MatRes['Similarity'] = KT
    if AddEdge:
        # MatRes['P1'] = G1.P
        # MatRes['P2'] = G2.P
        MatRes['Edges'] = G1.Edges
        MatRes['NEdges'] = G2.Edges
        MatRes['KQ'] = KQ
    MatRes['KP'] = KP
    MatRes['GT'] = range(NofNodes)    
    sio.savemat(TmpFName, MatRes)
    

    # Edges
    if AddEdge:
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

    # Triplet
    if AddTriplet:
        for ti in range(KT.shape[0]):
        #for ti in range(0):
            CTripletsVec = VecInt(3)
            CTripletsVec[0] = int(G1.Triplets[ti][0])
            CTripletsVec[1] = int(G1.Triplets[ti][1])
            CTripletsVec[2] = int(G1.Triplets[ti][2])
            CurrentNNZV = doubleArray(KT.shape[1])
            for xijk in range(KT.shape[1]):
                CurrentNNZV[xijk] = KT[ti][xijk]
            G.AddGenericGenericSparseFactor(CTripletsVec, nnzTripIdx, CurrentNNZV)

    G.AddAuctionFactor()
    
    return G, TmpFName;


def array_hash(d):
    """ dict => string hash """
    return '{' + ','.join(sorted(str(k) +':'+ str(d[k]) for k in range(len(d)))) + '}'


#~
# Constructing the higher order graphical model for super nodes produced by local modes search procedure.
#~
def ConstructSuperGraph(NofNodes, SPNodes, LargePositiveNumber = 50):
    NofStates = intArray(NofNodes)
    NofFactors = 0;
    
    for i in range(NofNodes):
        NofStates[i] = NofNodes
    G = CFactorGraph(NofNodes, NofStates)

    for (idx, SPNode) in SPNodes.iteritems():
        AllEntries = [e for e in SPNode.iterkeys()]
        e1 = eval(AllEntries[0]);
        clu = [n for n in e1.iterkeys()]
        t = sorted(clu)
        if len(clu) != len(t):
            print("Error! Wrong format for super nodes! A node appears multiple times.")
            print("The super node is")
            print(clu)
        #if (array_hash(t) in IsAdded):
        #    print("The super node is already added! Skipping it.")
        #    continue
        #else:
        #    IsAdded[array_hash(t)] = 1

        if(len(clu) > 1):
            CluSize = len(clu)
            CluVec = VecInt(CluSize)
            for i in range(CluSize):
                CluVec[i] = clu[i]
            NNZs = len(AllEntries)
            NNZVArray = [-v for v in SPNode.itervalues()]
            print (NNZVArray)
            NNZv = doubleArray(NNZs)
            NNZIdxVec = VecVecInt(NNZs)
            for xijk in range(len(AllEntries)):
                cAssignVec = VecInt(CluSize)
                cAssignDict = eval(AllEntries[xijk])
                cAssignArray = [xi for xi in cAssignDict.itervalues()]
                for xi in range(len(cAssignArray)):
                    cAssignVec[xi] = cAssignArray[xi]
                NNZv[xijk] = NNZVArray[xijk] + LargePositiveNumber
                NNZIdxVec[xijk] = cAssignVec
            G.AddGenericGenericSparseFactor(CluVec, NNZIdxVec, NNZv)
            NofFactors += 1
    G.AddAuctionFactor()
    return G, NofFactors

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
