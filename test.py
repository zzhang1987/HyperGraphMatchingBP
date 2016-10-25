import FactorBP as FB
import numpy as np


np.random.seed(3297461)


NofNodes = 3;
NofStates = FB.intArray(NofNodes);
for i in range(NofNodes):
    NofStates[i] = NofNodes

NofStates[0] = 3
NofStates[1] = 2
NofStates[2] = 3

G = FB.CFactorGraph(NofNodes, NofStates)



for i in range(NofNodes):
    Potentials = FB.doubleArray(NofStates[i])
    for j in range(NofStates[i]):
        Potentials[j] = np.random.normal();

    G.AddNodeBelief(i, Potentials);


EdgePotentials = FB.doubleArray(3 * 2);

for i in range(6):
    EdgePotentials[i] = np.random.normal();


mi = FB.doubleArray(3)
mj = FB.doubleArray(2)
mi[0] = 0
mj[0] = 0
mi[1] = 0
mj[1] = 0
mi[2] = 0



G.AddEdge(0, 1, EdgePotentials);

G.PrintFactorInfo();
for iter in range(100):
    G.UpdateMessages();

G.PrintFactorInfo();


