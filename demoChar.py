import scipy.io as sio
# import cv2
import numpy as np
import numpy.matlib
from FactorBP import *
from FactorBP.FactorGraph import *
from ChineseChar import *
#from BaBSolver import *
from drawMatches import *
from FactorBP.FindNearSol  import FindMModes
import time
import matplotlib
import multiprocessing as mp
from itertools import product
from scipy.optimize import linear_sum_assignment

# matplotlib.use('Agg')
# import matplotlib.pyplot as plt

def RunModelMModes((Fname, G, NofNodes, delta, N)):
    start_time = time.time()
    MModes = FindMModes(NofNodes, G, delta, N)
    time_dur = time.time() - start_time
    Fname = '%s_ID%d_NOus%d.pkl' % (Fname, idx, NofOus)
    f = open(Fname, "w")
    pickle.dump(MModes, f)
    pickle.dump(gTruth, f)
    pickle.dump(NofOus, f)
    pickle.dump(time_dur, f)
    f.close()
ErrorRate = np.zeros(10)


def ComputeMModesChinese(idx, d1, d2):
    idx1base = (idx) * 10
    if(id1 >= id2):
        return
    data1 = sio.loadmat('./data_chrct/' + str(idx1base + d1 + 1) + '.mat')
    G1 = np.array(data1['G'])
    I1 = np.array(data1['I'])
    Pt1 = np.array(data1['Pt'])

    data2 = sio.loadmat('./data_chrct/' + str(idx1base + d2 + 1) + '.mat')

    G2 = np.array(data2['G'])
    I2 = np.array(data2['I'])
    Pt2 = np.array(data2['Pt'])

    G = ConstructSparseG(G1, Pt1, G2, Pt2)
    gTruth = range(Pt1.shape[0])
    start_time = time.time()
    MModes = FindMModes(Pt1.shape[0], G, delta, N)
    time_dur = time.time() - start_time
    Fname = '%s_%d_%d_%d_char.pkl' % ('Char',idx, d1, d2)

    f = open(Fname, "w")
    pickle.dump(MModes, f)
    pickle.dump(gTruth, f)
    pickle.dump(NofOus, f)
    pickle.dump(time_dur, f)
    f.close()


def worker(f, inQ, outQ, kw=None):
    while True:
        x = inQ.get()
        if x is None:
            break
        outQ.put(f(x, **kw))

def pmap(f, tasks, n_jobs=mp.cpu_count(), **kw):
    n_jobs = min(n_jobs, len(tasks)) # XXX explicitly avoid zombies
    inQ, outQ = mp.Queue(), mp.Queue()
    procs = [mp.Process(target=worker, args=(f, inQ, outQ, kw)) for _ in range(n_jobs)]
    for p in procs:
        p.daemon = True
        p.start()
    for x in tasks:
        inQ.put(x)
    for _ in range(n_jobs):
        inQ.put(None)
    results = [outQ.get() for _ in range(len(tasks))]
    for p in procs:
        p.join()
    return results

Idxs = range(4)

for idx in range(4):
    cnt = 0;
    SumErrorRate = 0.0;
    idx1base = (idx ) * 10;
    AllTime = 0.0;
    delta = 8
    N = 300
    
    for d1 in range(10):
        for d2 in range(d1 + 1, 10):
            data1 = sio.loadmat('./data_chrct/' + str(idx1base + d1 + 1) + '.mat');
            cnt += 1
            G1 = np.array(data1['G'])
            I1 = np.array(data1['I'])
            Pt1 = np.array(data1['Pt'])

            data2 = sio.loadmat('./data_chrct/' + str(idx1base + d2 + 1) + '.mat')

            G2 = np.array(data2['G'])
            I2 = np.array(data2['I'])
            Pt2 = np.array(data2['Pt'])

            #G = ConstructSparseG(G1, Pt1, G2, Pt2)
            #start_time = time.time()
            #MModes = FindMModes(Pt1.shape[0], G, delta, N)
            #time_dur = time.time() - start_time
            G, K = ConstructSparseG(G1, Pt1, G2, Pt2)

            # G.SetVerbost(True)\
            #G = ConstructDenseG(G1, Pt1, G2, Pt2)
            # G.SetVerbost(True)
            G.SetVerbose(False)
            NofNodes = Pt1.shape[0]
            xhat = None
            bestv = 0
            K = 0.5 * (K + K.transpose())
            lastdual = 1e20
            for iter in range(1000):
                G.UpdateMessages()
                xc = G.GetCurrentDecode()
                xcv = G.CurrentPrimal()
                xc_mat = np.zeros([NofNodes, NofNodes])
                for xi in range(NofNodes):
                    xc_mat[xi][xc[xi]] = 1
                xc_vec = xc_mat.reshape([NofNodes * NofNodes, 1])
                last_cv = 0
                for iter in range(100):
                    c1 = K.dot(xc_vec.copy())
                    c1mat = c1.reshape([NofNodes, NofNodes])
                    row_inds, x1 = linear_sum_assignment(-c1mat)
                    x1_mat = np.zeros([NofNodes, NofNodes])
                    for xi in range(NofNodes):
                        x1_mat[xi][x1[xi]] = 1
                    x1_vec = x1_mat.reshape([NofNodes * NofNodes, 1])
                    cv = x1_vec.transpose().dot(c1)
                    if(abs(cv - last_cv) < 1e-6):
                        break
                    last_cv = cv
                    rv = x1_vec.transpose().dot(K.dot(x1_vec))
                    if(rv > bestv):
                        bestv = rv
                        xhat = x1_mat.argmax(axis=0)
                    xc_vec = x1_mat.reshape([NofNodes * NofNodes, 1])
                dual = G.DualValue()
                primal = bestv
                if(abs(dual - primal) < 1e-2):
                    break
                if(abs(dual - lastdual) < 1e-2):
                    break
                lastdual = dual

            #time_start = time.time()
            #G.Solve(1000)
            #time_dur = time.time() - time_start
            #res = BaBSolver(G,2,1000,0.0005, False);
            #print("Time=%.4f" % res.Time)
            #dG.SetVerbose(True)
            res = xhat
            corrected = np.ones([len(res), 1]);

            #AllTime += time_dur #res.Time

            decode = res
            ErrAssign = 0.0;
            #GTDecode = intArray(len(G1))
            for xi in range(len(G1)):
                if (decode[xi] != xi):
                    ErrAssign += 1
                    corrected[xi] = 0
                #GTDecode[xi] = xi
            #drawMatches(255 - I1,255 - I2, Pt1, Pt2, G1, G2, res.Decode, corrected)


            ErrorRate = ErrAssign / len(G1);
            SumErrorRate += ErrorRate
            #rValue = G.ComputeObj(GTDecode)

            # print(G.GetDecode())

    print("Char", idx, " Accuracy ", 1 - SumErrorRate / cnt) #, "Time ", AllTime/cnt)
    #plt.show()
