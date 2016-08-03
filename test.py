from FactorGraph import *;
import numpy as np;


np.random.seed(3297461);


NofNodes = 3;
NofStates = intArray(NofNodes);
for i in range(NofNodes):
    NofStates[i] = NofNodes

G = CFactorGraph(NofNodes, NofStates)


Potentials  = doubleArray(3);

for i in range(NofNodes):
    for j in range(NofNodes):
        Potentials[j] = np.random.normal();

    G.AddNodeBelief(i, Potentials);


EdgePotentials = doubleArray(NofNodes * NofNodes);

for i in range(NofNodes * NofNodes):
    EdgePotentials[i] = np.random.normal();


SparseEdgePotentials = doubleArray(NofNodes * NofNodes);

for i in range(NofNodes * NofNodes):
    SparseEdgePotentials[i] = 0;


    
NNZs = 3;

nnzIdx = intArray(2 * NofNodes);

nnzIdx[0] = 0
nnzIdx[1] = 0

nnzIdx[2] = 1
nnzIdx[3] = 2

nnzIdx[4] = 0
nnzIdx[5] = 1


for i in range(NNZs):
    idx = nnzIdx[2 * i] * NofNodes + nnzIdx[2 * i + 1]
    SparseEdgePotentials[idx] = np.random.normal();

mi = doubleArray(NofNodes);
mj = doubleArray(NofNodes);

for i in range(NofNodes):
    mi[i] = 0
    mj[i] = 0

G.AddEdge(0, 1, EdgePotentials);
G.AddSparseEdge(1, 1, SparseEdgePotentials, mi, mj, NNZs, nnzIdx)


G.PrintFactorInfo();

G.UpdateMessages();

G.PrintFactorInfo();
