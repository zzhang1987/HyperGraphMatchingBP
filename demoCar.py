import scipy.io as sio
# import cv2
import numpy as np
import numpy.matlib
from FactorBP import *
#from FactorGraph import *
from ConstructGM import *
#from BaBSolver import *
from sys import argv

def RunCar(id, outs, IsSparse):
    data = sio.loadmat('./Car/CarID' + str(int(id)) + '_nOut' + str(int(outs)) + '.mat');

    Edges1 = data['Edge1'];
    Edges2 = data['Edge2'];
    KP = data['KP']
    KQ = data['KQ']

    G = ConstructG(Edges1, Edges2, KP, KQ, IsSparse);
    G.SetVerbose(True)
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
    for nOus in range(1):
        for idx in range(1,2):
            print("Compute for Instance %d Outliers %d" % (idx, nOus))
            ResFname = './Car/CarIDRes' + str(idx) + '_nOut' + str(nOus) + '_' +str(True) + '.mat';
            res1 = RunCar(idx, nOus, True)
            Res = {}
            Res['res'] = res1
            print(res1)
            sio.savemat(ResFname, Res)
