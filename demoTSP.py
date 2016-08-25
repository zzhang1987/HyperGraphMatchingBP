import numpy as np
from FactorGraph import *
from FactorGraph import *
from ConstructGM import *
from BaBSolver import *
#from sys import argv

np.random.seed(123456)
size_pr = 20
c = np.random.uniform(0, 10.0, [size_pr, size_pr])
c = c + c.transpose();
c = c ;
NofNodes = int(size_pr)
NofStates = intArray(NofNodes)



oldc = c;
c = - c;

for i in range(size_pr):
    c[i][i] = -1000;
    NofStates[i] = NofNodes

EPotentials = doubleArray(NofNodes * NofNodes)

for i in range(NofNodes):
    for j in range(NofNodes):
        EPotentials[i * NofNodes + j] = c[i][j]

G = CFactorGraph(NofNodes, NofStates)
for i in range(size_pr):
    G.AddEdge(i, (i + 1)%size_pr, EPotentials)

G.AddAuctionFactor()
#G.SetVerbose(True)
#G.Solve(1000)
G.SetVerbose(False)
res = BaBSolver(G, 1000, 5, 0.0000005, True);
visited = np.zeros(size_pr)
IsValid = True
i = 0
cnt = 0
while(cnt < size_pr):
    i = res.Decode[i]
    if(visited[i]):
        IsValid = False
        break
    else:
        visited[i] = True
        cnt = cnt + 1

print(IsValid)
print(res.Decode)

