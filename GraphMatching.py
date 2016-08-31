from FactorGraph import *
from ConstructGM import *
from BaBSolver import *


class ModelGraph:
    def __init__(self, Img, Pt, Edge, LocalFeatures):
        """Construct a model graph
        :param Img the image
        :param Pt the keypoints
        :param Edge the edges between points
        :param LocalFeatures the local features
        """
        self.Img = Img
        self.Pt = Pt
        self.Edge = Edge
        self.LocalFeatures = LocalFeatures

    def ComputeKP(self, G2):
        """Compute node similarities
        :param G2 another model graph
        """
        raise NotImplementedError

    def ComputeKQ(self, G2):
        """Compute node similarities
        :param G2 another model graph
        """
        raise NotImplementedError

