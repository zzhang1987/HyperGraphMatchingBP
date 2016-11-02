import time

import numpy as np
class PQueueNode:
    Next = None;
    data = None;
    def __init__(self, data):
        self.data = data
        self.Next = None

class PQueue:
    Head = None
    def __init__(self):
        self.Head = None
    def insert(self, data):
        if self.Head is None:
            self.Head = PQueueNode(data)
            return
        nNode = PQueueNode(data)
        if (data < self.Head.data):
            nNode.Next = self.Head
            self.Head = nNode
            return
        p = self.Head

        while(p.Next != None):
            if data > p.data and p.Next.data > data:
                nNode.Next = p.Next
                p.Next = nNode
                return
            p=p.Next
        p.Next = nNode
    def extract_min(self):
        res = self.Head
        self.Head = self.Head.Next
        res.Next = None
        return res
    def isempty(self):
        return (self.Head == None)

class BaBNode:
    UB = np.nan
    LB = np.nan
    DualStore = None
    Node = None
    AssignMent = None
    Gap = None;
    def GetUB(self):
        return self.UB
    def __init__(self, DualStore, UB, LB, Node, AssignMent, Gap):
        self.DualStore = DualStore
        self.AssignMent = AssignMent
        self.Node = Node
        self.UB = UB
        self.LB = LB
        self.Gap = Gap

    def __del__(self):
        #print("Node deleted!");
        if self.DualStore is not None:
            del self.DualStore
            self.DualStore = None
        if self.Node is not None:
            del self.Node
        if self.AssignMent is not None:
            del self.AssignMent

    def __lt__(self, other):
        if(self.UB > other.UB):
            return True
        elif( (self.UB == other.UB) and (self.Gap < other.Gap)):
            return True
        return False

    def __gt__(self, other):
        if (self.UB <= other.UB):
            return True
        return False

class BaBRes:
    __slots__ = ["Value","Decode", "Time"]
    def __init__(self,value, decode, time):
        self.Value = value
        self.Decode = decode
        self.Time = time

def BaBSolver(G,outIter, inIter, maxGap, verbose):
    start_time = time.time()
    G.Solve(400)
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
        #for N in range(G.NofNodes()):
        #    Node = BaBNode(GStore,GB,LB,MFNode.Nodes, N + 1, MFNode.gap)
        #    Q.insert(Node)
        Node2 = BaBNode(GStore,GB,LB,MFNode.Nodes, -(MFNode.States + 1), MFNode.gap)
        Node1 = BaBNode(GStore,GB,LB,MFNode.Nodes, MFNode.States + 1, MFNode.gap)

        Q.insert(Node2)
        Q.insert(Node1)
        del cNode
        iter += 1

    time_dur = time.time() - start_time

    return BaBRes(GLB, xhat,time_dur)
