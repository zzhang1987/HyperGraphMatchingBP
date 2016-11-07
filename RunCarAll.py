import matlab.engine
eng = matlab.engine.start_matlab()
#import matlab.engine
import numpy as np;
import FactorBP as FB
import drawMatches as dm
import RunAlgorithm as RA
from Utils import LoadCar
import cPickle as pickle


def ComputeAccuracyPas(decode, gTruth, NofInliers ):
    Ccnt = 0
    for i in range(len(gTruth)):
        if((decode[i] == gTruth[i]) and (gTruth[i] < NofInliers)):
            Ccnt += 1
    return 1.0 * Ccnt / NofInliers
# eng = matlab.engine.start_matlab()
CarData = LoadCar()


AlgorithmNames=['Ours', 'OursBca', 'BCA', 'BCA-MP', 'BCA-IPFP', 'HGM', 'RRWHM', 'TM', 'OursPW', 'FGM']

SecondOrderMethods = ('Ours', 'FGM')
ThirdOrderMethods = ('Ours', 'OursBca', 'BCA', 'BCA-MP', 'BCA-IPFP', 'HGM', 'RRWHM', 'TM')

MaxNofOus = 20
NofInstances = 30
reload(RA)
AllAcc = dict()
AllRtime = dict()
AllObj = dict()

for NofOus in range(0,MaxNofOus+1):
    Accuracy = dict()
    Rtime = dict()
    Obj = dict()
    AllAcc[NofOus] = dict()
    AllRtime[NofOus] = dict()
    AllObj[NofOus] = dict()
    for idx in range(1, NofInstances + 1):
        car1 = CarData[idx]
        LocalFeature1 = car1['features1']
        LocalFeature2 = car1['features2']
        
        PT1 = LocalFeature1[:, 0:2]
        PT2 = LocalFeature2[:, 0:2]
        
        
        orientation1 = LocalFeature1[:, 8]
        orientation2 = LocalFeature2[:, 8]
        
        GT = car1['gTruth'][0]
        
        NofInliers = len(GT)
        CMaxNofOus = np.min([LocalFeature1.shape[0], LocalFeature2.shape[0]]) - NofInliers
        CNofOus = NofOus
        if(CNofOus > CMaxNofOus):
            CNofOus = CMaxNofOus
        NofNodes = CNofOus + NofInliers
        gTruth = np.random.permutation(NofNodes)
        PT1 = PT1[gTruth, :]
        orientation1 = orientation1[gTruth]
        MG1 = FB.MatchingGraph(PT1[0:NofNodes], orientation1[0:NofNodes])
        MG2 = FB.MatchingGraph(PT2[0:NofNodes], orientation2[0:NofNodes])
        for Type in ('pas', 'pasDisOnly'):
            for WithEdge in (True,False):
                if(WithEdge == True):
                    for methods in SecondOrderMethods:
                        FullMethods = (methods + 'WithEdge' + str(WithEdge) +
                                       'WithTriplet' + str(False) + Type)
                        print('Run %s' % FullMethods)
                        if(idx == 1):
                            Accuracy[FullMethods] = dict()
                            Rtime[FullMethods] = dict()
                            Obj[FullMethods] = dict()
                        decode,rtime,obj = RA.RunAlgorithm(MG1, MG2, WithEdge,
                                                           False,  Type, methods,
                                                           eng)
                        Accuracy[FullMethods][idx] = ComputeAccuracyPas(decode, gTruth, NofInliers)
                        Rtime[FullMethods][idx] = rtime
                        Obj[FullMethods][idx] = obj

                        Fname = ('Res/Car%d_Nous%d_' + FullMethods + '.pdf') % (idx, NofOus)
                        dm.drawMatchesWithOutlier(car1['I1'],car1['I2'],PT1[0:NofNodes],PT2[0:NofNodes],decode, gTruth, NofInliers, Fname)


                for methods in ThirdOrderMethods:
                
                    FullMethods = (methods + 'WithEdge' + str(WithEdge) +
                                   'WithTriplet' + str(True) + Type)
                    print('Run %s' % FullMethods)

                    if(idx == 1):
                        Accuracy[FullMethods] = dict()
                        Rtime[FullMethods] = dict()
                        Obj[FullMethods] = dict()
                    decode, rtime, obj = RA.RunAlgorithm(MG1, MG2, WithEdge,
                                                         True, Type, methods,
                                                         eng)
                    Accuracy[FullMethods][idx] = ComputeAccuracyPas(decode,
                                                                    gTruth,
                                                                    NofInliers)
                    Rtime[FullMethods][idx] = rtime
                    Obj[FullMethods][idx] = obj
                    Fname = ('Res/Car%d_Nous%d_' + FullMethods
                             + '.pdf') % (idx, NofOus)
                    dm.drawMatchesWithOutlier(car1['I1'], car1['I2'],
                                              PT1[0:NofNodes], PT2[0:NofNodes],
                                              decode, gTruth,
                                              NofInliers, Fname)
        
    AllAcc[NofOus] = Accuracy
    AllRtime[NofOus] = Rtime
    AllObj[NofOus] = Obj
                
f = open('CarRes.pkl', "w")
pickle.dump(AllAcc, f)
pickle.dump(AllRtime, f)
pickle.dump(AllObj, f)
f.close()

