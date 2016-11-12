import FactorBP as FB
import matplotlib
matplotlib.use('Agg')
import skimage
import skimage.io as imageio
from skimage.feature import (match_descriptors, corner_harris,
                             corner_peaks, ORB, plot_matches)
import matplotlib.pyplot as plt

I1 = imageio.imread('FourthOrder/boat/img1.pgm')
imageio.imshow(I1)
I2 = imageio.imread('FourthOrder/boat/img2.pgm')
imageio.imshow(I2)
descriptor_extractor = ORB(n_keypoints=80)

descriptor_extractor.detect_and_extract(I1)
keypoints1 = descriptor_extractor.keypoints
descriptors1 = descriptor_extractor.descriptors

descriptor_extractor.detect_and_extract(I2)
keypoints2 = descriptor_extractor.keypoints
descriptors2 = descriptor_extractor.descriptors

matches12 = match_descriptors(descriptors1, descriptors2, cross_check=True)

fig = plt.figure()
ax = plt.axes(aspect=1)
plot_matches(ax, I1, I2, keypoints1, keypoints2, matches12)

MG1 = FB.MatchingGraph(keypoints1, descriptors1)
MG2 = FB.MatchingGraph(keypoints2, descriptors2) 

G,MFname = FB.ConstructMatchingModel(MG1, MG2, 'syn', True, False)
G.solve(1000)
p=range(80)
nmatches = np.array([np.array(p), G.GetDecode()]).transpose()
