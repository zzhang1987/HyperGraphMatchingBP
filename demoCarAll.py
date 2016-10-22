import matlab.engine
import numpy as np;
from sklearn.neighbors import KDTree
import FactorBP as FB
from scipy.spatial import Delaunay
import scipy.io as sio
from Utils import *


def ComputeAccuracyPas(decode, NofInliers):
    Ccnt = 0
    for i in xrange(NofInliers):
        if (decode[i] == i):
            Ccnt += 1
    return 1.0 * Ccnt / NofInliers


# eng = matlab.engine.start_matlab()
CarData = LoadCar()
Accurancy = np.zeros(30)
AccurancyPW = np.zeros(30)

NofOus = 0
for idx in xrange(1, 31):
    car1 = CarData[idx]
    LocalFeature1 = car1['features1']
    LocalFeature2 = car1['features2']
    PT1 = LocalFeature1[:, 0:2]
    PT2 = LocalFeature2[:, 0:2]
    orientation1 = LocalFeature1[:, 8]
    orientation2 = LocalFeature2[:, 8]
    GT = car1['gTruth'][0]
    NofInliers = len(GT)
    NofNodes = NofOus + NofInliers
    MG1 = FB.MatchingGraph(PT1[0:NofNodes], orientation1[0:NofNodes])
    MG2 = FB.MatchingGraph(PT2[0:NofNodes], orientation2[0:NofNodes])

    G1 = FB.ConstructMatchingModel(MG1, MG2, 'pas', True)
    G2 = FB.ConstructMatchingModel(MG1, MG2, 'pas', False)
    G1.SetVerbose(False)
    G2.SetVerbose(False)
    res = FB.BaBSolver(G1, 600, 5, 0.005, False)
    resPW = FB.BaBSolver(G2, 600, 5, 0.005, False)
    Accurancy[idx - 1] = ComputeAccuracyPas(res.Decode, NofInliers)

    print('%d %f %f' % (idx, ComputeAccuracyPas(res.Decode, NofInliers), ComputeAccuracyPas(resPW.Decode, NofInliers)))
    # print(res.Time)
    # print(res.Value)

    AccurancyPW[idx - 1] = ComputeAccuracyPas(resPW.Decode, NofInliers)

    # print(ComputeAccuracyPas(resPW.Decode, NofInliers))
    # print(resPW.Time)
    # print(resPW.Value)
print(Accurancy.mean())
print(AccurancyPW.mean())