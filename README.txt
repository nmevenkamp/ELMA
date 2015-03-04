### ELMA - ELectron Micrograph Analysis

Version 1.0

ELMA is a small tool for denoising and analyzing STEM images. The main application are low-dose high-resolution STEM images of inorganic materials that are extremely corrupted by Poisson noise. Currently, ELMA is focused on periodic crystals, but we intend to extend it to more complex structures in the future.

Executables are available for Windows and MacOSX at http://nmevenkamp.github.io/ELMA.

ELMA supports TIFF (integer & float) input/output compatible with DigitalMicrograph. The list of supported formats is:
* .tif, .tiff, .dm3 (read only), .dm4 (read only), .png, .pgm, .q2bz

Images must be in units of electron counts / pixel.

LICENSE: ELMA is distributed under the terms of the Common Development and Distribution License (see LICENSE.txt).  
Particularly, you can download and use this software free of cost.  
   
We appreciate any feedback on your experience with ELMA. In case you encounter any problems when using this software, please don't hesitate to contact us: mevenkamp@aices.rwth-aachen.de

If you publish results generated using ELMA, please refer to our related work:

[1] N. Mevenkamp, P. Binev, W. Dahmen, P.M. Voyles, A.B. Yankovich, and B. Berkels: Poisson Noise Removal from High-Resolution STEM Images based on periodic Block-matching, Advanced Structural and Chemical Imaging, 2015.

[2] N. Mevenkamp, A.B. Yankovich, P.M. Voyles, and B. Berkels: Non-local Means
    for Scanning Transmission Electron Microscopy Images and Poisson Noise based
    on Adaptive Periodic Similarity Search and Patch Regularization, in
    J. Bender, A. Kuijper, T. von Landesberger, H. Theisel, and P. Urban (Eds.),
    VMV 2014: Vision, Modeling & Visualization, pp. 63–70.
    Eurographics Association, Darmstadt, Germany (2014). 