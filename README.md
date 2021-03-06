### ELMA - ELectron Micrograph Analysis

__Version 1.0__

ELMA is a small tool for **denoising** and analyzing STEM images. The main application are **low-dose high-resolution STEM images** of inorganic materials that are extremely corrupted by **Poisson noise**. Currently, ELMA is focused on periodic crystals, but we intend to extend it to more complex structures in the future.
___
Executables are available for **Windows** and **MacOSX** at [http://nmevenkamp.github.io/ELMA](http://nmevenkamp.github.io/ELMA).

ELMA supports **TIFF** (integer & float) **input/output compatible with DigitalMicrograph**. The list of supported formats is:
* **.tif**, **.tiff**, **.dm3** (read only), **.dm4** (read only), .png, .pgm, .q2bz

Images must be in units of **electron counts / pixel**.
___
**LICENSE**: ELMA is distributed under the terms of the [Common Development and Distribution License](LICENSE.txt).  
Particularly, you can **download and use** this software **free of cost**.  
   
We appreciate any feedback on your experience with ELMA. In case you encounter any problems when using this software, please don't hesitate to contact us: [mevenkamp@aices.rwth-aachen.de](mailto:mevenkamp@aices.rwth-aachen.de)
___
If you publish results generated using ELMA, please refer to our related work:

**[1]** N. Mevenkamp, P. Binev, W. Dahmen, P.M. Voyles, A.B. Yankovich, and B. Berkels:
    **Poisson noise removal from high-resolution STEM images based on periodic block matching**.
    Advanced Structural and Chemical Imaging 2015, 1:3. [ [PDF](http://www.ascimaging.com/content/pdf/s40679-015-0004-8.pdf) | [DOI](http://dx.doi.org/10.1186/s40679-015-0004-8) | [BIB](http://www.aices.rwth-aachen.de:8080/~mevenkamp/papers/BIRS2015/MeBiDaVoYaBe2015.bib) ]

**[2]** N. Mevenkamp, A.B. Yankovich, P.M. Voyles, and B. Berkels: **Non-local Means**
    **for Scanning Transmission Electron Microscopy Images and Poisson Noise based**
    **on Adaptive Periodic Similarity Search and Patch Regularization**, in
    J. Bender, A. Kuijper, T. von Landesberger, H. Theisel, and P. Urban (Eds.),
    VMV 2014: Vision, Modeling & Visualization, pp. 63–70.
    Eurographics Association, Darmstadt, Germany, 2014. [ [PDF](http://www.aices.rwth-aachen.de:8080/~mevenkamp/papers/VMV2014/MeYaVoBe2014.pdf) | [DOI](http://dx.doi.org/10.2312/vmv.20141277) | [BIB](http://www.aices.rwth-aachen.de:8080/~mevenkamp/papers/VMV2014/MeYaVoBe2014.bib) ]