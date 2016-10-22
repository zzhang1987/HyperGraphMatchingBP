import numpy as np
from scipy.spatial import Delaunay
import scipy.io as sio


def computeFeatureSimple(Points, T):
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
        F[idx] = vecX[((idx + 1) % 3)] * vecY[idx] - vecY[((idx + 1) % 3)] * vecX[idx]

    return F


def CreateTensorHouseDelaunay(P1, P2, bpermute):
    res = {}
    NP1 = P1.shape[0]
    NP2 = P2.shape[0]
    if(bpermute):
        res['GT'] = np.random.permutation(P2.shape[0])
    else:
        res['GT'] = np.array(range(P2.shape[0]))
    P1 = P1[res['GT'], :]
    tri1 = Delaunay(P1)
    tri2 = Delaunay(P2)
    t1 = tri1.simplices
    t2 = PermunateTriplets(tri2.simplices)
    
    #Because of the super symmetric, we only need to permunate t2
    
    Feature1 = computeFeatureSingle(P1,t1)
    Feature2 = computeFeatureSingle(P2,t2)
    
    distMat = ComputeFeatureDistance(Feature1, Feature2)
    
    dist = np.exp(- (distMat / np.mean(distMat) ))
    
    res['Triplets'] = t1
    res['NTriplets'] = t2
    res['Similarity'] = dist
    return res

def IndicesToVec(indices, NofNodes, NofStates):
    res = VecInt(NofNodes)
    res[2] = indices % NofStates
    indices /= NofStates
    res[1] = indices % NofStates
    indices /= NofStates
    res[0] = indices
    return res

def ComputeAccuracy(X1, GT):
    Ccnt = 0
    for b in range(GT.shape[0]):
        if(int(X1[b]) == int(GT[b])):
            Ccnt += 1
    res = Ccnt * 1.0 / GT.shape[0]
    return res



def computeFeatureSingle(P, T):
    Feature = np.zeros([T.shape[0], 3])
    for i in range(T.shape[0]):
        Feature[i]  = computeFeatureSimple(P, T[i])
    return Feature
def ComputeFeatureDistance(F1,F2):
    res = np.zeros([F1.shape[0], F2.shape[0]])
    for i in range(F1.shape[0]):
        for j in range(F2.shape[0]):
            res[i][j] = np.linalg.norm(F1[i]-F2[j])
        
    return res

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

def LoadCar():
    BaseSuffix = './Cars_and_Motorbikes_Graph_Matching_Datasets_and_Code/Data_for_Cars_and_Motorbikes/Data_Pairs_Cars/pair_%d.mat'
    NofPairs = 30
    res = {}
    for i in range(1, NofPairs+1):
        res[i] = sio.loadmat(BaseSuffix % i)

    return res

def LoadMotor():
    BaseSuffix = './Cars_and_Motorbikes_Graph_Matching_Datasets_and_Code/Data_for_Cars_and_Motorbikes/Data_Pairs_Motorbikes/pair_%d.mat'
    NofPairs = 20
    res = {}
    for i in range(1, NofPairs+1):
        res[i] = sio.loadmat(BaseSuffix % i)

    return res
