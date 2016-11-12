import numpy as np
import FactorGraph as FG
import BaBSolver as BSolver
from scipy.spatial.distance import hamming
import cPickle as pickle
import MatchingGraph as MG
import FactorGraph as FG


def SequentialDivMBest(NofNodes, G, delta, N, MaxIter = 1000):
    res = dict()
    G.Solve(1000)
    DualStore = G.StoreDual()
    G.ResetMax()
    Xarray = FG.intArray(NofNodes)
    for i in range(N):
        G.Solve(1000)
        G.ResetMax()
        X = G.GetDecode()
        for j in range(NofNodes):
            Xarray[j] = int(X[j])
        X1 = np.array(X, dtype=np.int32)
        T = G.StoreDual()
        G.ReStoreDual(DualStore)
        v = G.ComputeObj(Xarray)
        res[X1.tostring()] = v
        G.ReStoreDual(T)
        for i in range(NofNodes):
            G.AddValue(i, int(X[i]), -delta)
    return res

def ParallelDivMBest(NofNodes, G,  NofSolutions, MaxIter = 1000):
    res = dict()



    for iter in range(MaxIter):
        G.UpdateMessages()
        dual = G.DualValue()
        print("Iter = %d, Dual = %f " % (iter, dual))
    FinalDecode = G.GetDecode()
    for i in range(NofSolutions):
        res[i] = np.zeros(NofNodes, dtype=np.int32)
        for xi in range(NofNodes):
            res[i][xi] = FinalDecode[i][xi]
    return res



def RunDataPDiverse((Fname, data, idx, NofOus, NofSolutions, delta)):
    car1 = data[idx]
    N = NofSolutions
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
    if (CNofOus > CMaxNofOus):
        CNofOus = CMaxNofOus
    NofNodes = CNofOus + NofInliers
    gTruth = np.random.permutation(NofNodes)
    PT1 = PT1[gTruth, :]

    orientation1 = orientation1[gTruth]
    MG1 = MG.MatchingGraph(PT1[0:NofNodes], orientation1[0:NofNodes])
    MG2 = MG.MatchingGraph(PT2[0:NofNodes], orientation2[0:NofNodes])
    G, MFname = MG.ConstructMatchingModelPDiverse(MG1, MG2,
                                       'pas', False, True, N, delta)

    res = ParallelDivMBest(NofNodes, G, N)
    Fname = '%s_ID%d_NOus%d_Delta_%f_PDiverse.pkl' % (Fname, idx, NofOus, delta)

    f = open(Fname, "w")
    pickle.dump(res, f)
    pickle.dump(gTruth, f)
    pickle.dump(NofOus, f)
    f.close()


def RunDataSDiverse((Fname, data, idx, NofOus, delta)):
    car1 = data[idx]
    N = 10
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
    
    MDiverse = SequentialDivMBest(NofNodes, G, delta, N)
    
    Fname = '%s_ID%d_NOus%d_Delta_%f_SDiverse.pkl' % (Fname, idx, NofOus, delta)

    f = open(Fname, "w")
    pickle.dump(MDiverse, f)
    pickle.dump(gTruth, f)
    pickle.dump(NofOus, f)
    f.close()

