import scipy.io as sio
# import cv2
import numpy as np
import numpy.matlib
from FactorGraph import *
from ChineseChar import *

ErrorRate = np.zeros(10)

for idx in range(1):
    cnt = 0;
    SumErrorRate = 0.0;
    idx1base = (idx) * 10;

    for d1 in range(1,2):
        for d2 in range(d1 + 1, 3):

            data1 = sio.loadmat('./data_chrct/' + str(idx1base + d1 + 1) + '.mat');
            cnt += 1
            G1 = np.array(data1['G'])
            I1 = np.array(data1['I'])
            Pt1 = np.array(data1['Pt'])

            data2 = sio.loadmat('./data_chrct/' + str(idx1base + d2 + 1) + '.mat');

            G2 = np.array(data2['G'])
            I2 = np.array(data2['I'])
            Pt2 = np.array(data2['Pt'])

            G = ConstructSparseG(G1, Pt1, G2, Pt2)
            # G.SetVerbost(True)\
            #dG = ConstructDenseG(G1, Pt1, G2, Pt2)
            # G.SetVerbost(True)
            G.SetVerbose(True)
            #dG.SetVerbose(True)
            G.Solve(100)
            #dG.Solve(100)
            decode = G.GetDecode()
            ErrAssign = 0.0;
            #GTDecode = intArray(len(G1))
            for xi in range(len(G1)):
                if (decode[xi] != xi):
                    ErrAssign += 1;
                #GTDecode[xi] = xi

            ErrorRate = ErrAssign / len(G1);
            SumErrorRate += ErrorRate
            #rValue = G.ComputeObj(GTDecode)

            # print(G.GetDecode())
    print("Char", idx, " Error Rate ", 1 - SumErrorRate / cnt)
