import numpy as np;
import FactorBP as FB
from Utils import LoadCar
from FactorBP.FindNearSol import FindModes
import ipdb


def ComputeAccuracyPas(decode, gTruth, NofInliers ):
    Ccnt = 0
    for i in range(NofInliers):
        if(decode[i] == gTruth[i]):
            Ccnt += 1
    return 1.0 * Ccnt / NofInliers


ipdb.set_trace()
CarData = LoadCar()

NofOus = 10
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

G3,MFname = FB.ConstructMatchingModel(MG1, MG2, 'pas', AddEdge = True, AddTriplet = False)

#XMAP = FB.BaBSolver

#X0 = np.random.permutation(NofNodes)
#delta = 6
#X = FindModes(NofNodes, G3, X0, delta)

