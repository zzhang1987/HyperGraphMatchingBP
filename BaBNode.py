import numpy as np

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
        if (self.UB < other.UB):
            return True
        elif ((self.UB == other.UB) & (self.Gap >= other.Gap)):
            return True
        return False