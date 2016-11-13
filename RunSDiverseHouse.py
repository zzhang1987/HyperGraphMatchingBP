from Utils import LoadHouse
from FactorBP.DiverseMBest import RunDataSDiverseHouse
import numpy as np

import multiprocessing as mp
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
    n_jobs = min(n_jobs, len(tasks))  # XXX explicitly avoid zombies
    inQ, outQ = mp.Queue(), mp.Queue()
    procs = [mp.Process(target=worker, args=(f, inQ, outQ, kw))
             for _ in range(n_jobs)]
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


HouseData, HouseImage = LoadHouse()

ImageIs = range(0,111)
BaseLines = range(10, 101, 10)
Lambdas = np.arange(0.1,4,0.1)

Nous = 5

pmap(RunDataSDiverseHouse, [('House', HouseData,
                             ImageI, baseline, Nous, 10, delta) for
                            (delta, baseline, ImageI) in product(Lambdas,
                                                                 BaseLines,
                                                                 ImageIs)],
     n_jobs=32)

