import numpy as np
import FactorGraph as FG
import BaBSolver as BSolver
from scipy.spatial.distance import hamming
import cPickle as pickle
import MatchingGraph as MG
import FactorGraph as FG
import HouseUtils as HU
import time


def ComputeAccuracyPas(decode, gTruth, NofInliers):
    Ccnt = 0
    for i in range(len(gTruth)):
        if((decode[i] == gTruth[i]) and (gTruth[i] < NofInliers)):
            Ccnt += 1
    return 1.0 * Ccnt / NofInliers


def SequentialDivMBest(NofNodes, G, delta, N,
                       MaxIter=1000, gTruth=None, NofInliers=None):
    res = dict()
    G.Solve(1000)
    DualStore = G.StoreDual()
    G.ResetMax()
    Xarray = FG.intArray(NofNodes)
    CurrentAccuracy = 0
    CV = 0
    for i in range(N):
        G.ResetMax()
        G.Solve(1000)
        X = G.GetDecode()
        for j in range(NofNodes):
            Xarray[j] = int(X[j])
        X1 = np.array(X, dtype=np.int32)
        T = G.StoreDual()
        G.ReStoreDual(DualStore)
        v = G.ComputeObj(Xarray)
        if(gTruth is not None and NofInliers is not None):
            CA = ComputeAccuracyPas(X, gTruth, NofInliers)
            if(CA > CurrentAccuracy):
                CurrentAccuracy = CA
                CV = v
            print('Iter = %d, Accuracy = %f, CA = %f, CV = %f' % (i,
                                                                  CurrentAccuracy, CA, v))
        
        res[X1.tostring()] = v
        G.ReStoreDual(T)
        for i in range(NofNodes):
            G.AddValue(i, int(X[i]), -delta)
    return res

def ParallelDivMBest(NofNodes, G,  NofSolutions, MaxIter = 10000):
    res = dict()
    for iter in range(MaxIter):
        G.UpdateMessages()
        dual = G.DualValue()
        #print("Iter = %d, Dual = %f " % (iter, dual))
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
    start_time = time.time()
    res = ParallelDivMBest(NofNodes, G, N)
    time_dur = time.time() - start_time

    Fname = '%s_ID%d_NOus%d_Delta_%f_PDiverse.pkl' % (Fname, idx, NofOus, delta)

    f = open(Fname, "w")
    pickle.dump(res, f)
    pickle.dump(gTruth, f)
    pickle.dump(NofOus, f)
    pickle.dump(time_dur, f)
    f.close()

    
def RunDataPDiverseHouse((Fname, HouseData, ImageI,
                          baseline, NofOus, NofSols, delta)):
    MG1, MG2, gTruth = HU.GenerateDataHouse(HouseData, ImageI,
                                            baseline, NofOus)
    if MG1 is None:
        return None

    G, MFname = MG.ConstructMatchingModelPDiverse(MG1, MG2,
                                                  'pas', False,
                                                  True, NofSols, delta)
    start_time = time.time()
    res = ParallelDivMBest(30 - NofOus, G, NofSols)
    time_dur = time.time() - start_time

    Fname = '%s_ID_%d_BaseLine%d_NOus%d_Delta_%f_PDiverse.pkl' % (Fname, ImageI, baseline, NofOus, delta)

    f = open(Fname, "w")
    pickle.dump(res, f)
    pickle.dump(gTruth, f)
    pickle.dump(NofOus, f)
    pickle.dump(time_dur, f)
    f.close()

    
def RunDataSDiverseHouse((Fname, HouseData, ImageI,
                          baseline, NofOus, NofSols, delta)):
    
    MG1, MG2, gTruth = HU.GenerateDataHouse(HouseData, ImageI,
                                            baseline, NofOus)
    if MG1 is None:
        return None
    G, MFname = MG.ConstructMatchingModel(MG1, MG2,
                                          'cmu', True, False)
    start_time = time.time()
    MDiverse = SequentialDivMBest(30 - NofOus, G, delta,
                                  NofSols)
    time_dur = time.time() - start_time

    Fname = '%s_ID_%d_BaseLine%d_NOus%d_Delta_%f_SDiverse.pkl' % (Fname, ImageI, baseline, NofOus, delta)

    f = open(Fname, "w")
    pickle.dump(MDiverse, f)
    pickle.dump(gTruth, f)
    pickle.dump(NofOus, f)
    pickle.dump(time_dur, f)
    f.close()

    
def RunDataSDiverse((Fname, data, idx, NofOus, delta)):
    car1 = data[idx]
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
    G, MFname = MG.ConstructMatchingModel(MG1, MG2, 'pas', False, True)
    start_time = time.time()
    MDiverse = SequentialDivMBest(NofNodes, G, delta, N, gTruth=gTruth,
                                  NofInliers=NofInliers)
    time_dur = time.time() - start_time
    
    Fname = '%s_ID%d_NOus%d_Delta_%f_SDiverse.pkl' % (Fname, idx, NofOus, delta)

    f = open(Fname, "w")
    pickle.dump(MDiverse, f)
    pickle.dump(gTruth, f)
    pickle.dump(NofOus, f)
    pickle.dump(time_dur, f)
    f.close()

