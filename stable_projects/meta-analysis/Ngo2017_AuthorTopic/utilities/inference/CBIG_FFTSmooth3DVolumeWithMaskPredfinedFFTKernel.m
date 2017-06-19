function smoothedV = FFTSmooth3DVolumeWithMaskPredfinedFFTKernel(V, local_mask, fftK)
% smoothedV = FFTSmooth3DVolumeWithMaskPredfinedFFTKernel(V, local_mask, fftK)

% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

tmpV = V;
tmpV(local_mask == 0) = 0;


smoothedV = CBIG_FFTSmooth3DWithPredefinedKernel(tmpV, fftK);  
smoothedV(local_mask ~= 0) = bsxfun(@times, smoothedV(local_mask ~= 0), 1./local_mask(local_mask~=0));
smoothedV(local_mask == 0) = V(local_mask == 0);
