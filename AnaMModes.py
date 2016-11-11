from itertools import product
import sys
import numpy as np
import scipy.optimize as so
import scipy.stats as ss
from scipy.spatial.distance import hamming

if(sys.version_info >= (3,0)):
    import _pickle as pickle
else:
    import cPickle as pickle


def ComputeAccuracyPas(decode, gTruth, NofInliers):
    Ccnt = 0
    for i in range(len(gTruth)):
        if((decode[i] == gTruth[i]) and (gTruth[i] < NofInliers)):
            Ccnt += 1
    return 1.0 * Ccnt / NofInliers


Index = range(1,31)
NOus = range(0,21)
Fname = 'Car'
MeanAccMap = np.zeros(len(NOus))
MeanAccMax = np.zeros(len(NOus))
for NofOus in NOus:
    AccuracyMax = np.zeros(len(Index))
    AccuracyMap = np.zeros(len(Index))
    
    for idx in Index:
        FileName = 'Diverse/%s_ID%d_NOus%d_Delta_0.600000_SDiverse.pkl' % (Fname, idx, NofOus)
        f = open(FileName, "rb")
        Modes = pickle.load(f)
        gTruth = pickle.load(f)
        NofOus = pickle.load(f)
        f.close()
        
        NofNodes = gTruth.shape[0]
        AllV = np.array(list(Modes.values()))
        MV = AllV.mean()
        NofInliers = NofNodes - NofOus
        Prob = np.zeros([NofNodes, NofNodes])
        MAP = -100
        MAPAcc = 0
        MAPLabel = 0
        MaxAcc = 0
        MaxAccV = 0
        for key in Modes.keys():
            cLabel = np.fromstring(key, dtype=np.int32)
            cValue = Modes[key]
            if(cValue > MAP):
                MAP = cValue
                MAPAcc = ComputeAccuracyPas(cLabel, gTruth, NofInliers)
                MAPLabel = cLabel
        for key in Modes.keys():
            cLabel = np.fromstring(key, dtype=np.int32)
            cValue = Modes[key]

            if(ComputeAccuracyPas(cLabel, gTruth, NofInliers) > MaxAcc):
                MaxAcc = ComputeAccuracyPas(cLabel, gTruth, NofInliers)
                MaxAccV = cValue
        #print('MAP Accuracy: %f' % MAPAcc)
        #print('Max Accuracy: %f' % MaxAcc)
        AccuracyMap[idx - 1] = MAPAcc
        AccuracyMax[idx - 1] = MaxAcc
    MeanAccMap[NofOus] = AccuracyMap.mean()
    MeanAccMax[NofOus] = AccuracyMax.mean()

    print('MeanAccMAP: %f MeanAccMax: %f' % (MeanAccMap[NofOus], MeanAccMax[NofOus]))
