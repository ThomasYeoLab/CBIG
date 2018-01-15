% compile all mex cpp scripts 
% Code has been updated (4/4/2011) to work with 64 bit compilers
mex -O -largeArrayDims GibbsSamplerLDA.cpp
mex -O -largeArrayDims GibbsSamplerAT.cpp
mex -O -largeArrayDims GibbsSamplerHMMLDA.cpp
mex -O -largeArrayDims GibbsSamplerLDACOL.cpp
mex -O -largeArrayDims GibbsSamplerLDA_NEWDOCS.cpp
mex -O -largeArrayDims binarysearchstrings.c

