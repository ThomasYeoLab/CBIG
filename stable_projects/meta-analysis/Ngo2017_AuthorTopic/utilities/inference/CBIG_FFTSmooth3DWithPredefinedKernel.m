function smoothedV = FFTSmooth3DWithPredefinedKernel(V, fftK)
  % smoothedV = FFTSmooth3DWithPredefinedKernel(V, fftK)

  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  fftV = fftn(fftshift(V)); 
  tmp = fftV .* fftK;
  smoothedV = ifftn(tmp);
