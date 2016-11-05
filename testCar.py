import numpy as np;
import FactorBP as FB
from Utils import LoadCar
from FactorBP.FindNearSol import FindModes

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
car1 = CarData[idx]
LocalFeature1 = car1['features1']
LocalFeature2 = car1['features2']

PT1 = LocalFeature1[:, 0:2]
PT2 = LocalFeature2[:, 0:2]



orientation1 = LocalFeature1[:, 8].copy()
orientation2 = LocalFeature2[:, 8].copy()

GT = car1['gTruth'][0] 

NofInliers = len(GT)
NofNodes = NofOus + NofInliers
gTruth = np.random.permutation(NofNodes)
PT1 = PT1[gTruth, :]
norientation = orientation1[gTruth].copy()
orientation1[0:NofNodes] = norientation


PT1cp = PT1.copy()
PT2cp = PT2[0:NofNodes].copy()



MG1 = FB.MatchingGraph(PT1[0:NofNodes], norientation)
MG2 = FB.MatchingGraph(PT2[0:NofNodes], orientation2[0:NofNodes])

#G2 = FB.ConstructMatchingModel(MG1, MG2, 'pas', True)
#G2.SetMinDualDecrease(1e-6)
#G2.Solve(200)
#res =  FB.BaBSolver(G2, 100, 5, 0.005, True)
#print(ComputeAccuracyPas(res.Decode, gTruth, NofInliers))
print(gTruth)
#print(G2.GetDecode())

G3,MFname = FB.ConstructMatchingModel(MG1, MG2, 'pasDis', AddEdge = True, AddTriplet = False)

#XMAP = FB.BaBSolver

#X0 = np.random.permutation(NofNodes)
#delta = 6
#X = FindModes(NofNodes, G3, X0, delta)
DualStore = G3.StoreDual()

res = FB.BaBSolver(G3, 100, 10, 0.005, True)


delta1 = 12
delta2 = 8


for i in range(100):
    X1 = GenRandomLabel(res.Decode, NofNodes)
    G3.ReStoreDual(DualStore)
    G3.ResetMax()
    X2 = FindModes(NofNodes, G3, X1, delta2)

    print(X2)

    print(ComputeAccuracyPas(X2, gTruth, NofInliers))

