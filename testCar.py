import numpy as np;
import FactorBP as FB
from Utils import LoadCar, LoadHouse
from FactorBP.FindNearSol import RunDataMModes, RunDataMModesHouse
from FactorBP.DiverseMBest import RunDataSDiverse, RunDataPDiverse, RunDataSDiverseHouse, RunDataPDiverseHouse


def ComputeAccuracyPas(decode, gTruth, NofInliers ):
    Ccnt = 0
    for i in range(NofInliers):
        if(decode[i] == gTruth[i]):
            Ccnt += 1
    return 1.0 * Ccnt / NofInliers
def GenRandomLabel(X0, delta1):
    t1 = np.random.permutation(X0.shape[0])
    t2 = t1[0:delta1]
    t3 = np.random.permutation(t2)
    X1 = X0.copy()
    X1[t2] = X0[t3]
    return X1

CarData = LoadCar()

NofOus = 5
idx = 1
np.random.seed(123456)
HouseData,HouseImg = LoadHouse()

RunDataSDiverseHouse(('House',HouseData,0,10,5,10,0.1))
RunDataPDiverseHouse(('House',HouseData,0,10,5,10,0.1))
RunDataMModesHouse(('House', HouseData, 0, 10, 5))

#RunDataPDiverse(('Car', CarData, 1, 0, 10, 2.5))

#RunDataSDiverse(('Car', CarData, 1, 0, 0.5))

#RunDataMModes(('Car', CarData, 1, 0))


