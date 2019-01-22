function smoothedV = ...
    CBIG_AuthorTopic_FFTSmooth3DVolumeWithMaskPredfinedFFTKernel(V, ...
    mask, fftK)
% smoothedV = ...
%   CBIG_AuthorTopic_FFTSmooth3DWithPredefinedKernel(V, ...
%   mask, fftK)
%
% Smooth a 3D volume using Fast Fourier Transform (FFT) with a predefined
% smoothing kernel and a mask
%
% Input:
%  - V: the original 3D volume.
%  - mask: mask with the same size as V. All voxels with value 0 in the
%     in the mask are set to 0 in the smoothed volume. All voxels with
%     value > 0 in the mask are not affected.
%  - fftK: the smoothing kernel.
% Output:
%  - smoothedV: the output volume after smoothing and masking.
%
% Example:
%   smoothedV = ...
%     CBIG_AuthorTopic_FFTSmooth3DVolumeWithMaskPredfinedFFTKernel(V, ...
%     mask, presetKernel)
%   Smooth the 3D volume V using the kernel presetKernel and return the
%   final smoothed volume as smoothedV. All voxels with value 0 in mask
%   are set to 0 in smoothedV.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

tmpV = V;
tmpV(mask == 0) = 0;


smoothedV = CBIG_AuthorTopic_FFTSmooth3DVolumeWithPredefinedKernel(tmpV, fftK);
smoothedV(mask ~= 0) = bsxfun(@times, smoothedV(mask ~= 0), 1./mask(mask~=0));
smoothedV(mask == 0) = V(mask == 0);
