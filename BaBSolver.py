import numpy as np
import copy
from FactorGraph import *
from BaBNode import *
from PQueue import *

import time
from FibonacciHeap import FibHeap

class BaBRes:
    __slots__ = ["Value","Decode", "Time"]
    def __init__(self,value, decode, time):
        self.Value = value
        self.Decode = decode
        self.Time = time

def BaBSolver(G,outIter, inIter, maxGap, verbose):
    start_time = time.time()
    time_dur = None
    #G.Solve(20)
    Q = PQueue()
    iter = 1
    #xhat = G.GetDecode()
    GUB = 1e20;#G.DualValue()
    GLB = -1e20#G.PrimalValue()
    root = BaBNode(None,GUB,GLB,None,None,None)

    Q.insert(root)
    while(Q.isempty() == False):
        cNode = Q.extract_min()
        if (iter >= outIter):
            break
        if (cNode == None):
            print("No more node in Queue, Exact solution found.\n")
            break
        cNode = cNode.data
        if(verbose):
            print("BaBIter=%d, GUB = %f, GLB = %f, Gap = %.2f, Time=%.4f" % (iter, GUB, GLB, np.fabs(GUB - GLB)/np.fabs(GLB), time.time() - start_time ))

        if(GUB < GLB):
            del cNode
            break;
        if (np.fabs(GUB - GLB) / np.fabs(GLB) < maxGap):
            del cNode
            break;
        if(cNode.DualStore != None):
            G.ReStoreDual(cNode.DualStore)
        if(cNode.Node != None):
            G.SetDecode(cNode.Node, cNode.AssignMent)
        GUB = cNode.UB
        G.Solve(inIter)
        LB = G.PrimalValue()
        GB = G.DualValue()

        if(LB > GLB):
            GLB = LB;
            xhat = G.GetDecode()
        if (np.fabs(LB - GB) < 1e-5):
            iter = iter + 1
            del cNode
            continue
        if (GB < GLB):
            iter = iter + 1
            del cNode
            continue

        MFNode = G.FindMostFracNodes()
        GStore = G.StoreDual()
        Node1 = BaBNode(GStore,GB,LB,MFNode.Nodes, MFNode.States + 1, MFNode.gap)
        Node2 = BaBNode(GStore,GB,LB,MFNode.Nodes, -(MFNode.States + 1), MFNode.gap)

        Q.insert(Node1)
        Q.insert(Node2)
        del cNode
        iter += 1

    time_dur = time.time() - start_time

    return BaBRes(GLB, xhat,time_dur)
