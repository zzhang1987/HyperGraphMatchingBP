import numpy as np;
import FactorBP as FB
from Utils import LoadCar
from FactorBP.FindNearSol import RunDataMModes
from FactorBP.DiverseMBest import RunDataSDiverse, RunDataPDiverse


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
RunDataPDiverse(('Car', CarData, 1, 0, 10, 0.5))

RunDataSDiverse(('Car', CarData, 1, 0, 0.5))

RunDataMModes(('Car', CarData, 1, 0))


