#import matlab.engine
import numpy as np;
from sklearn.neighbors import KDTree
import FactorBP as FB
from scipy.spatial import Delaunay
import scipy.io as sio
from Utils import *

# eng = matlab.engine.start_matlab()
CarData = LoadCar()

car1 = CarData[5]
LocalFeature1 = car1['features1']
LocalFeature2 = car1['features2']
PT1 = LocalFeature1[:, 0:2]
PT2 = LocalFeature2[:, 0:2]
orientation1 = LocalFeature1[:, 8]
orientation2 = LocalFeature2[:, 8]
GT = car1['gTruth'][0]
NofInliers = len(GT)
NofOus = 0
NofNodes = NofInliers + NofOus

MG1 = FB.MatchingGraph(PT1[0:NofNodes], orientation1[0:NofNodes])
MG2 = FB.MatchingGraph(PT2[0:NofNodes], orientation2[0:NofNodes])

G1 = FB.ConstructMatchingModel(MG1, MG2, 'pas', True)
G2 = FB.ConstructMatchingModel(MG1, MG2, 'pas', False)
G1.SetVerbose(True)
G2.SetVerbose(False)
G1.Solve(1000)
res = FB.BaBSolver(G1, 600, 10, 0.005, True)
resPW = FB.BaBSolver(G2, 600, 10, 0.005, True)
print(res.Decode)
print(res.Time)
print(res.Value)

print(resPW.Decode)
print(resPW.Time)
print(resPW.Value)