import scipy.io as sio
import cv2
data1 = sio.loadmat('./Ch/data_chrct/1.mat');

G1 = data1['G']
I1 = test['I']
Pt1 = test['Pt']

data2 = sio.loadmat('./Ch/data_chrct/2.mat');

G2 = data2['G']
I2 = data2['I']
Pt1 = data2['Pt']
