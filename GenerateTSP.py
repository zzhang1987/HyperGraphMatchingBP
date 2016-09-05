import numpy as np

def GenerateEuclideanTSP(n):
    distMat = np.zeros([n,n]);
    for i in range(n):
        distMat[i,i] = 10000;
    points = np.random.uniform(0.0, 10.0, [n, 2]);

    for i in range(n):
        for j in range(i+1, n):
            distMat[i,j] = np.linalg.norm(points[i] - points[j])
            distMat[j,i] = distMat[i,j]

    return distMat

def GenerateSymTSP(n):
    distMat = np.random.uniform(0.0, 10.0, [n, n]);
    distMat = 0.5 * (distMat.transpose() + distMat)
    
    for i in range(n):
        distMat[i,i] = 10000;

    return distMat

def GenerateASymTSP(n):
    distMat = np.random.uniform(0.0, 10.0, [n, n]);
    
    for i in range(n):
        distMat[i,i] = 10000;

    return distMat

    
