function smoothedPhi = CBIG_AuthorTopic_SmoothPhi(params)
% smoothedPhi = CBIG_AuthorTopic_SmoothPhi(params)
%
% Spatially smooth the parameters of the variational distribution phi.
%
% Input:
%  - params: struct containing parameters of the author-topic model and
%            the Collapsed Variational Bayes (CVB) algorithm.
%            The spatial smoothing is performed with the mask saved at
%            params.maskPath. The 3D spherical Gaussian smoothing kernel is
%            defined with the radius params.smoothKernelSize, and standard
%            deviation parms.smoothKernelSigma set up by
%            CBIG_AuthorTopic_SetupParameters.m
% Output:
%  - smoothedPhi: spatially smoothed paramaeters phi of the variational
%                 distribution.
%
% Example:
%   params = CBIG_AuthorTopic_SetupParameters(23, 2, '/AT_outputs', ...
%              '/data/setGeneratedThoughtData.mat', 50)
%   Set up the CVB algorithm's parameters. The CVB algorithm is randomly
%   initialized with seed 23. The author-topic model is assumed to have 2
%   components. The CVB algorithm has its output saved at /AT_outputs.
%   The algorithm uses input from /data/selfGeneratedThoughtData.mat.
%   The CVB algorithm runs for 50 iterations before checking for
%   convergence.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  smoothedPhi = params.phi;
  brainMask   = MRIread(params.maskPath);
  nonzero_vol  = zeros(1, params.V, 'single');
  masked_vol   = zeros(1, params.V, 'single');
  % filter

  tmpVol = zeros(size(brainMask.vol, 1), size(brainMask.vol, 2), size(brainMask.vol, 3), 'single');
  mask    = zeros(size(brainMask.vol, 1), size(brainMask.vol, 2), size(brainMask.vol, 3), 'single');

  kernel = CBIG_AuthorTopic_3DGaussianLPF(params.smoothKernelSize, params.smoothKernelSigma);
  paddedKernel = CBIG_AuthorTopic_Pad3DSmoothingKernel(kernel, size(tmpVol));
  fftK = fftn(paddedKernel);
  for e = 1:params.E
    fociIndices              = params.fociIndices{e};
    masked_vol(~fociIndices) = 0;
    masked_vol(fociIndices)  = 1;
    mask(brainMask.vol(:) ~= 1) = 0;
    mask(brainMask.vol(:) == 1) = masked_vol;

    smoothedMask = CBIG_AuthorTopic_FFTSmooth3DVolumeWithPredefinedKernel(mask, fftK);
    smoothedMask(mask ~= 1) = 0;
    for k = 1:params.K
      for t = 1:params.expByNumTasks(e)
        nonzero_vol(~fociIndices) = 0;
        nonzero_vol(fociIndices)  = squeeze(params.phi{e}(:, k, t));
        tmpVol(brainMask.vol(:) ~= 1) = 0;
        tmpVol(brainMask.vol(:) == 1) = nonzero_vol;

        tmpVol = CBIG_AuthorTopic_FFTSmooth3DVolumeWithMaskPredfinedFFTKernel(tmpVol, smoothedMask, fftK);

        tmp = tmpVol(brainMask.vol(:) == 1);
        smoothedPhi{e}(:, k, t) = tmp(fociIndices);
      end;
    end;
    denom = sum(sum(smoothedPhi{e},3),2);
    smoothedPhi{e} = bsxfun(@times, smoothedPhi{e}, 1./denom);
    smoothedPhi{e}(denom == 0,:,:) = 0;
  end;



function paddedKernel = CBIG_AuthorTopic_Pad3DSmoothingKernel(kernel, newSize)
% paddedKernel = CBIG_AuthorTopic_Pad3DSmoothingKernel(kernel, newSize)
%
% Pad a 3D kernel with zeros.
%
% Input:
%  - kernel : original 3D kernel. Each dimension of the kernel is
%             preferably an odd integer.
%  - newSize: a 3x1 vector specifying the kernel's size after padding.
%             Each dimension of the padded kernel is preferably an odd
%             integer.
% Output:
%  - paddedKernel: the new 3D kernel after padding.
%
% Example:
%   smallKernel = rand(3 3 3);
%   paddedKernel = CBIG_AuthorTopic_Pad3DSmoothingKernel(smallKernel, ...
%     [11 11 11]);
%   Pad the 3x3x3 kernel smallKernel of the size with zeros in order
%     to get a new kernel of the size 11x11x11
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  paddedKernel = zeros(newSize(1), newSize(2), newSize(3));
  kSize = size(kernel);
  paddedKernel(...
      uint64(newSize(1)/2-kSize(1)/2+1):uint64(newSize(2)/2+kSize(1)/2), ...
      uint64(newSize(2)/2-kSize(2)/2+1):uint64(newSize(2)/2+kSize(2)/2), ...
      uint64(newSize(3)/2-kSize(3)/2+1):uint64(newSize(3)/2+kSize(3)/2)) = kernel;



function h = CBIG_AuthorTopic_3DGaussianLPF(dim, sigma)
% h = CBIG_AuthorTopic_3DGaussianLPF(dim, sigma)
%
% Generate a rotationally symmetric 3D Gaussian lowpass filter with
% standard deviation sigma (in pixels) and dimension specified by dim
% (in pixels).
%
% Input:
%  - dim: a scalar or a 1x3 vector specifying
%         the dimensions of the kernel. Default value: [3 3 3]
%  - sigma: standard devivation of the 3D Gaussian. Default value: 0.65
% Output:
%  - h: a 3D Gaussian lowpass filter.
%
% Example:
%   kernel = CBIG_AuthorTopic_3DGaussianLPF(31, 2)
%   Return a 31x31x31 pixels rotationally symmetric 3D Gaussian low-pass
%   filter with standard deviation of 2 pixels.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  if nargin>0,
    if ~(all(size(dim)==[1 1]) || all(size(dim)==[1 3])),
       error(id('InvalidFirstInput'),'The first parameter must be a scalar or a 1-by-3 size vector.');
    end
    if length(dim)==1, siz = [dim dim dim]; else siz = dim; end
  end

  if nargin<1, siz = [3 3 3]; end
  if nargin<2, std = .65; else std = sigma; end
  [x,y,z] = meshgrid(-(siz(2)-1)/2:(siz(2)-1)/2, -(siz(1)-1)/2:(siz(1)-1)/2, -(siz(3)-1)/2:(siz(3)-1)/2);
  h = exp(-(x.*x + y.*y + z.*z)/(2*std*std));
  h = h/sum(h(:));
