import numpy as np;
from sklearn.neighbors import KDTree
import FactorBP as FB
from scipy.spatial import Delaunay
import scipy.io as sio
from Utils import *
import time
from IPython.display import clear_output
import drawMatches as dm
import cPickle as pickle


def ComputeAccuracyPas(decode, gTruth, NofInliers ):
    Ccnt = 0
    for i in range(len(gTruth)):
        if((decode[i] == gTruth[i]) and (gTruth[i] < NofInliers)):
            Ccnt += 1
    return 1.0 * Ccnt / NofInliers

def RunAlgorithm(MG1, MG2, WithEdge, WithTriplet, Type, AlgoName, eng):
    NofNodes = MG1.P.shape[0]
    G,MFname = FB.ConstructMatchingModel(MG1, MG2, Type, WithTriplet, WithEdge)
    Gvis,MFname1 = FB.ConstructMatchingModel(MG1, MG2, Type, WithTriplet, WithEdge)
    cDecode = FB.intArray(NofNodes)

    if(AlgoName == 'Ours'):
        res1 = FB.BaBSolver(G, 120, 10, 0.005, False)
        Decode = res1.Decode
        RTime = res1.Time
        Obj = res1.Value
        return Decode, RTime, Obj
    if(AlgoName == 'OursBca'):
        res1 = FB.BaBSolver(G, 120, 10, 0.005, False)
        start_time = time.time()
        ResForBca = sio.loadmat(MFname)
        X0 = np.zeros(NofNodes)
        X0Vec = res1.Decode
        for i in xrange(NofNodes):
            X0[i] = X0Vec[i]
        ResForBca['X0'] = X0
        sio.savemat(MFname, ResForBca)
        resOursBCA = eng.runBcagm(MFname, nargout=3)
        time_dur = time.time() - start_time
            
        for i in range(NofNodes):
            cDecode[i] = int(resOursBCA[1][0][i])
        if(res1.Value < Gvis.ComputeObj(cDecode)):
            Decode = np.array(resOursBCA[1][0])
            RTime = res1.Time + time_dur
            Obj = Gvis.ComputeObj(cDecode)
        else:
            Decode = res1.Decode
            RTime = res1.Time
            Obj = res1.Value
        return Decode,RTime,Obj
    res = None
    if(AlgoName == 'BCA'):
        res = eng.runBcagm(MFname,nargout=3)
    if(AlgoName == 'BCA-MP'):
        res = eng.runBcagmQuad1(MFname,1,nargout=3)
    if(AlgoName == 'BCA-IPFP'):
        res = eng.runBcagmQuad1(MFname,2,nargout=3)
    if(AlgoName == 'HGM'):
        res = eng.runHGM(MFname,nargout=3)
    if(AlgoName == 'RRWHM'):
        res = eng.runRRWHM(MFname,nargout=3)
    if(AlgoName == 'TM'):
        res = eng.runTensorMatching(MFname,nargout=3)
    if(AlgoName == 'FGM'):
        res = eng.runFGM(MFname, nargout=3)

    if(res is not None):
        Decode = np.array(res[1][0])
        RTime = res[0]
        for i in range(NofNodes):
            cDecode[i] = int(res[1][0][i])
        Obj = Gvis.ComputeObj(cDecode)
    else:
        Decode = None
        RTime = None
        Obj = None

    return Decode, RTime, Obj
            
        
            




