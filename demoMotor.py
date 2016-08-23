import scipy.io as sio
# import cv2
import numpy as np
import numpy.matlib
from FactorGraph import *
from ConstructGM import *
from BaBSolver import *
from sys import argv

def RunMotor(id, outs, IsSparse):
    data = sio.loadmat('./Motor/MotorID' + str(int(id)) + '_nOut' + str(int(outs)) + '.mat');

    Edges1 = data['Edge1'];
    Edges2 = data['Edge2'];
    KP = data['KP']
    KQ = data['KQ']

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
    file_name, idx, nOus, IsSparse = argv
    idx = int(argv[1])
    nOus = int(argv[2])
    IsSparse = bool(int(argv[3]))
    print("Compute for Instance %d Outliers %d " % (idx, nOus))
    print(IsSparse)
    ResFname = './Motor/MotorIDRes' + str(idx) + '_nOut' + str(nOus) +'_' +str(IsSparse) + '.mat';
    res1 = RunMotor(idx, nOus, IsSparse)
    Res = {}
    Res['res'] = res1
    sio.savemat(ResFname, Res)
