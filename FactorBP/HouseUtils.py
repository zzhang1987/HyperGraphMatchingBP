import numpy as np
from MatchingGraph import MatchingGraph


def GenerateDataHouse(HouseData, ImageI, baseline, NofOus):
    np.random.seed(ImageI * baseline)
    if(ImageI + baseline > 110):
        return None, None, None
    PT1 = np.copy(HouseData[ImageI])
    PT2 = np.copy(HouseData[ImageI+baseline])
    NofNodes = 30
    rNofNodes = NofNodes - NofOus

    NPT1 = PT1.copy()
    PT1[:, 0] = NPT1[:, 1]
    PT1[:, 1] = NPT1[:, 0]
    NPT2 = PT2.copy()
    PT2[:, 0] = NPT2[:, 1]
    PT2[:, 1] = NPT2[:, 0]
    
    gTruth = np.random.permutation(NofNodes)
    PT1 = PT1[gTruth, :]

    gTruth2 = np.random.permutation(NofNodes)
    gTruth2 = np.array(sorted(gTruth2[0:rNofNodes]))
        
    for i in range(rNofNodes):
        N1 = gTruth[i]
        if(N1 not in gTruth2):
            gTruth[i] = -1
        elif N1 in gTruth2:
            gTruth[i] = np.argwhere(gTruth2 == N1)
                
    nPT1 = PT1[0:rNofNodes]
    nPT2 = PT2[gTruth2]
    PF1 = np.zeros([NofNodes - NofOus, 1])
    PF2 = np.zeros([NofNodes - NofOus, 1])

    MG1 = MatchingGraph(nPT1[0:rNofNodes], PF1[0:rNofNodes])
    MG2 = MatchingGraph(nPT2[0:rNofNodes], PF2[0:rNofNodes])

    return MG1, MG2, gTruth
