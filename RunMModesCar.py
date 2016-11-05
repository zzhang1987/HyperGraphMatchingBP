import numpy as np;
from Utils import LoadCar
from FactorBP.FindNearSol import RunDataMModes

import multiprocessing as mp
from contextlib import contextmanager
from itertools import product


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



CarData = LoadCar()
Index = range(1,31)
NOus = range(0,21)

pmap(RunDataMModes, [('Car', CarData, idx, NofOus) for (idx,NofOus) in product(Index,NOus) ], n_jobs=32)
