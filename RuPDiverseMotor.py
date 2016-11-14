from Utils import LoadCar, LoadMotor
from FactorBP.FindNearSol import RunDataMModes
from FactorBP.DiverseMBest import RunDataSDiverse, RunDataPDiverse

import multiprocessing as mp
from itertools import product
import numpy as np

#~
# parallelism
#~

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



CarData = LoadMotor()
Index = range(1,21)
NOus = range(0,21)
Deltas = np.arange(0.1,1.01,0.1)
pmap(RunDataPDiverse, [('Motor', CarData, idx, NofOus,10, Delta) for (idx,NofOus,Delta) in product(Index,NOus,Deltas) ], n_jobs=32)
