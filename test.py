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

G.AddEdge(0, 1, EdgePotentials);

G.PrintFactorInfo();

G.UpdateMessages();

G.PrintFactorInfo();
