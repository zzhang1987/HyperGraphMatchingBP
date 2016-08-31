import scipy.io as sio
# import cv2
import numpy as np
import numpy.matlib
from FactorGraph import *
from ChineseChar import *
from BaBSolver import *
from drawMatches import *
import matplotlib.pyplot as plt


ErrorRate = np.zeros(10)

for idx in range(4):
    cnt = 0;
    SumErrorRate = 0.0;
    idx1base = (idx ) * 10;
    AllTime = 0.0;
    for d1 in range(10):
        for d2 in range(d1 + 1, 10):
            
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
            G.SetVerbose(False)
            #G.Solve(5)
            #GStore = G.StoreDual();
            res = BaBSolver(G,600,5,0.005, False);
            #print("Time=%.4f" % res.Time)
            #dG.SetVerbose(True)

            corrected = np.ones([len(res.Decode), 1]);

            AllTime += res.Time

            decode = res.Decode
            ErrAssign = 0.0;
            #GTDecode = intArray(len(G1))
            for xi in range(len(G1)):
                if (decode[xi] != xi):
                    ErrAssign += 1
                    corrected[xi] = 0
                #GTDecode[xi] = xi
            drawMatches(255 - I1,255 - I2, Pt1, Pt2, G1, G2, res.Decode, corrected)


            ErrorRate = ErrAssign / len(G1);
            SumErrorRate += ErrorRate
            #rValue = G.ComputeObj(GTDecode)

            # print(G.GetDecode())

    print("Char", idx, " Accuracy ", 1 - SumErrorRate / cnt, "Time ", AllTime/cnt)
    plt.show()