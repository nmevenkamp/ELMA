#!/bin/bash
export CC=/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang
export CXX=/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang
cmake -G"Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DUSE_BOOST=1 -DDYNAMIC_LINKING=0 -DPARSE_GCC_ERRORS=1 -DUSE_BLAS=1 -DUSE_LAPACK=1 -DUSE_QT=1 -DBUILD_AND_USE_KISSFFT=1 -DUSE_PNG=1 -DUSE_TIFF=1 -DUSE_C++11=1 -DELMA_DEPLOY=1 ../quocmesh