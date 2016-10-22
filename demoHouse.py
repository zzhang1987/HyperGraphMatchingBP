import numpy as np;
from sklearn.neighbors import KDTree 
from FactorBP import *
from FactorBP.FactorGraph import *
from scipy.spatial import Delaunay
import scipy.io as sio
import matlab.engine

def LoadHouse():
    res = np.zeros([111, 30, 2])
    for i in range(1,112):
        res[i-1] = np.loadtxt('data/cmum/house/label/house%d' %i)
    return res

def computeFeatureSimple(Points, T):
    vecX = np.zeros(3)
    vecY = np.zeros(3)
    F = np.zeros(3)
    if((T[0]==T[1]) or (T[0]==T[2]) or (T[1]==T[2])):
        F = -10 * np.ones(3)
        return F
    for idx in xrange(3):
        vecX[idx] = Points[T[(idx + 1)%3]][0] - Points[T[idx]][0]
        vecY[idx] = Points[T[(idx + 1)%3]][1] - Points[T[idx]][1]
        length = np.linalg.norm([vecX[idx], vecY[idx]])
        if(length != 0):
            vecX[idx] /= length
            vecY[idx] /= length
        else:
            vecX[idx] = 0
            vecY[idx] = 0
  
    
    for idx in xrange(3):
        F[idx] = vecX[((idx+1)%3)]*vecY[idx]-vecY[((idx+1)%3)]*vecX[idx]
    
    return F

def computeFeatureSingle(P, T):
    Feature = np.zeros([T.shape[0], 3])
    for i in xrange(T.shape[0]):
        Feature[i]  = computeFeatureSimple(P, T[i])
    return Feature

def computeFeature(P1, P2, T1 ):
    res = {}
    F1 = np.zeros([T1.shape[1], 3])
    NP2 = P2.shape[0]
    F2 = np.zeros([NP2 * NP2 * NP2, 3])
    for i in xrange(T1.shape[1]):
        F1[i] = computeFeatureSimple(P1, T1[:,i])
    Fcnt = 0
    for i in xrange(NP2):
        for j in xrange(NP2):
            for k in xrange(NP2):
                F2[Fcnt] = computeFeatureSimple(P2, [i,j,k])
                Fcnt = Fcnt + 1
    res['feat1'] = F1
    res['feat2'] = F2
    
    return res

def CreateTensorHouse(P1,P2,bpermute):
    res = {}
    NP1 = P1.shape[0]
    NP2 = P2.shape[0]
    if(bpermute):
        res['GT'] = np.random.permutation(P2.shape[0])
    else:
        res['GT'] = np.array(range(P2.shape[0]))
    P1 = P1[res['GT'], :]
    
    nT = NP1 * NP2
    
    t1 = np.floor(np.random.rand(3, nT) * NP1)
    while(True):
        probFound=False;
        for i in range(3):
            ind=(t1[i,:]==t1[(i+1)%3,:])
            if(np.sum(ind)!=0):
                idxs = np.nonzero(ind)
                t1[i][idxs]=np.floor(np.random.rand(1,len(idxs[0]))*NP1);
                probFound=True;
        if(probFound == False):
            break;
    
    tri1 = Delaunay(P1)
    tri2 = Delaunay(P2)
    t1 = tri1.simplices
    t2 = tri2.simplices
    Feature = computeFeature(P1, P2, t1)
    
    kdt = KDTree(Feature['feat2'], metric='euclidean')
    nNN = 2000;
    [dist, indices] = kdt.query(Feature['feat1'], k=nNN, return_distance=True)
    
    dist = np.exp(- (dist / np.mean(dist) ))
    
    res['Triplets'] = t1
    res['NTriplets'] = indices
    res['Similarity'] = dist
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

    
def ComputeFeatureDistance(F1,F2):
    res = np.zeros([F1.shape[0], F2.shape[0]])
    for i in xrange(F1.shape[0]):
        for j in xrange(F2.shape[0]):
            res[i][j] = np.linalg.norm(F1[i]-F2[j])
        
    return res

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
    for b in xrange(GT.shape[0]):
        if(int(X1[b]) == int(GT[b])):
            Ccnt += 1
    if(Ccnt == 1):
        print(X1, GT)
    res = Ccnt * 1.0 / GT.shape[0]
    #print(Ccnt)
    return res

np.random.seed(123456)
HouseData = LoadHouse()


baseline = 90
start = 0
end = 111
TensorBCAAccuracy = np.zeros(end - baseline )
Accuracy = np.zeros(end - baseline )
Rtime = np.zeros(end - baseline )
eng = matlab.engine.start_matlab()


for ImageI in range(start, end - baseline):
    res = CreateTensorHouseDelaunay(np.copy(HouseData[ImageI]), np.copy(HouseData[ImageI+baseline]), True)
    sio.savemat('Temp.mat', res)

    NofNodes = 30
    NofStates = intArray(30)
    for i in range(30):
        NofStates[i] = 30
    G = CFactorGraph(NofNodes, NofStates)
    for i in range(res['Triplets'].shape[0]):
        T = res['Triplets'][i]
        T1 = VecInt(3)
        T1[0] = int(T[0])
        T1[1] = int(T[1])
        T1[2] = int(T[2])
        NNZIs = res['NTriplets']
        NNZVs = res['Similarity'][i]
        NNZIVecs = VecVecInt(NNZIs.shape[0])
        NNZVVecs = doubleArray(NNZIs.shape[0])
        for xi in xrange(NNZIs.shape[0]):
            cNTriplets = VecInt(3)
            cNTriplets[0] = int(NNZIs[xi][0])
            cNTriplets[1] = int(NNZIs[xi][1])
            cNTriplets[2] = int(NNZIs[xi][2])
            NNZIVecs[xi] = cNTriplets
            NNZVVecs[xi] = NNZVs[xi]
        G.AddGenericGenericSparseFactor(T1, NNZIVecs, NNZVVecs)

    G.SetVerbose(False)
    #G.Solve(1000)
    G.AddAuctionFactor()
    G.Solve(1000)
    res1 = BaBSolver(G, 50, 5, 0.005, False)
    #print(res1.Decode)
    #print(res['GT'])
    Accuracy[ImageI] = ComputeAccuracy(res1.Decode, res['GT'])
    Rtime[ImageI] = res1.Time
    print('Running Time: %f Accurancy: %f Valus: %f' %(Rtime[ImageI], Accuracy[ImageI], res1.Value))
    resBag = eng.runBcagm(nargout=3)
    #print(resBag)
    TensorBCAAccuracy[ImageI] =  ComputeAccuracy(resBag[1][0], res['GT']);
    print('TensorBCA: %f' % ComputeAccuracy(resBag[1][0], res['GT']))
    #print(res)

    

    
#print(Accuracy)
print(np.mean(Accuracy))
print(np.mean(TensorBCAAccuracy))
#print(np.mean(Rtime))

