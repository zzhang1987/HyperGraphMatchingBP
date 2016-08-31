import scipy.io as sio
# import cv2
import numpy as np
import numpy.matlib
from FactorGraph import *
from ChineseChar import *
from BaBSolver import *
from MatlabWrapper import *
ErrorRate = np.zeros(10)

for idx in range(4):
    cnt = 0;
    SumErrorRate = 0.0;
    idx1base = (idx ) * 10;
    AllTime = 0.0;
    for d1 in range(10):
        for d2 in range(d1 + 1, 10):
            Fname='Char/Char' + str(idx+1) + '_' + str(d1+1) + '_' + str(d2+1) + '.mat'
            res = RunPBP(Fname, True)
            Res = {}
            Res['res'] = res
            SaveFname = 'Char/Char' + str(idx+1) + '_' + str(d1+1) + '_' + str(d2+1) + 'SpRes.mat'
            sio.savemat(SaveFname, Res)
            res = RunPBP(Fname, False)
            Res = {}
            Res['res'] = res
            SaveFname = 'Char/Char' + str(idx + 1) + '_' + str(d1 + 1) + '_' + str(d2 + 1) + 'DeRes.mat'
            sio.savemat(SaveFname, Res)

