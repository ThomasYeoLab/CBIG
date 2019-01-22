function smoothedV = ...
    CBIG_AuthorTopic_FFTSmooth3DVolumeWithPredefinedKernel( ...
    V, fftK)
% smoothedV = ...
%   CBIG_AuthorTopic_FFTSmooth3DVolumeWithPredefinedKernel( ...
%   V, fftK)
%
% Smooth a 3D volume using Fast Fourier Transform (FFT) with a predefined
% smoothing kernel
%
% Input:
%  - V: the original 3D volume
%  - fftK: the smoothing kernel
% Output:
%  - smoothedV: the output smoothed volume
%
% Example:
%   smoothedV = CBIG_AuthorTopic_FFTSmooth3DVolumeWithPredefinedKernel( ...
%     V, presetKernel)
%   Smooth the 3D volume V using the kernel presetKernel and return the
%   smoothed volume as smoothedV
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  fftV = fftn(fftshift(V));
  tmp = fftV .* fftK;
  smoothedV = ifftn(tmp);
