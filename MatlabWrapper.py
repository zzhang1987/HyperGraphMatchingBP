import scipy.io as sio
# import cv2
import numpy as np
import numpy.matlib
from FactorGraph import *
from ConstructGM import *
from BaBSolver import *
from sys import argv

def RunPBP(Fname, IsSparse):
    data = sio.loadmat(Fname);

    Edges1 = data['Edge1'];
    Edges2 = data['Edge2'];
    KP = data['KP']
    KQ = data['KQ']

    #print(KP)
    
    G = ConstructG(Edges1, Edges2, KP, KQ, IsSparse);
    G.SetVerbose(False)
    # G.Solve(5)
    # GStore = G.StoreDual();
    res = BaBSolver(G, 600, 5, 0.005, False);
    res1 = []
    res1.append(res.Decode)
    res1.append(res.Value)
    res1.append(res.Time)
    del G
    return res1
if __name__ == '__main__':
    res = RunPBP(True)
    print(res)
