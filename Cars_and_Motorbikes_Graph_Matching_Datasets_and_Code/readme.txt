
Author: Marius Leordeanu
Date: May 24, 2013
Contact: leordeanu@gmail.com    


Cars and Motorbikes Graph Matching Datasets


The current archive contains the data (image pairs, correspondences and features) for the Cars and Motorbikes datasets introduced and used on graph matching in the following two papers:

[1]. Marius Leordeanu, Martial Hebert and Rahul Sukthankar, “Integer Projected Fixed Point Method for Graph Matching and MAP Inference”, NIPS 2009

and 

[2]. Marius Leordeanu, Rahul Sukthankar and Martial Hebert, Unsupervised Learning for Graph Matching, IJCV, 2012


This data is for research purposes only and, if used for publication, the above mentioned papers should be cited. 

_________________________________________________________________________________________

A note on the eventual experimental comparison: 

Paper [1] mentioned above introduces IPFP, an efficient method for graph matching 
which works both as a stand-alone-method (starting from a uniform solution), 
as well as a discretization procedure for obtaining a discrete/hard solution 
given the continuous output of any other graph matching algorithm (used for intialization).

When using these datasets please also compare against IPFP - uniform (starting from uniform solution) 
as well as IPFP  - SM (starting from the continuous solution of spectral matching [3]). 
Matlab code of these methods can be found here:  www.imar.ro/clvp/code.php . 
They are also included in the current archive in  ./Code_for_Spectral_Matching_and_IPFP.

[3]. Marius Leordeanu and Martial Hebert, “A Spectral Technique for Correspondence Problems using Pairwise Constraints”, ICCV 2005

________________________________________________________________________________________

Information on the datasets:

There are two image pair sequences, one for cars (30 image pairs) and another
one for motorbikes (20 image pairs). Here is a short description
to help you use it:

The data contains everything, from images, to contours and features that were used.

Please load the data and you will find the features in features1 (for
image 1) and features2 (for image 2).

They contain both inliers (features that have a correspondence, chosen manually)
and outliers (chosen randomly, for which there is probably no correct
correspondence).

A few details about how to interpret the data:

features1(i,1:2) - location of feature i from image 1.
features1(i,9) - orientation of the contour normal for feature i from image 1

These are the only pieces the information used in experiments from [1] and [2].

The ground truth can be found in the vector gTruth, it goes from 1 to
N, meaning that there are N inliers and for each inlier features1(i,:)
from image 1, corresponds to features2(i,:) from image 2.

Basically the first N features are ordered such that it is easy to
find the correspondences between images. The rest of the features, in
both features1 and features2 from N+1 to the end there are outliers.
If you visualize the locations you will probably see better what I
mean.

If you find any problems or bugs with the code or data then please let me know,
i would greatly appreciate it ! :) My email is: leordeanu@gmail.com
