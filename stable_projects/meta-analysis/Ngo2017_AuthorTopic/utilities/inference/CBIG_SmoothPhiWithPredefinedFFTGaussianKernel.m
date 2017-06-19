function smoothed_phi = CBIG_SmoothPhiWithPredefinedFFTGaussianKernel(params)
  % smoothed_phi = CBIG_SmoothPhiWithPredefinedFFTGaussianKernel(params)

  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  smoothed_phi = params.phi;
  brain_mask   = MRIread(params.mask_path);
  nonzero_vol  = zeros(1, params.V, 'single');
  masked_vol   = zeros(1, params.V, 'single');
  % filter
  
  tmp_vol = zeros(size(brain_mask.vol, 1), size(brain_mask.vol, 2), size(brain_mask.vol, 3), 'single');
  mask    = zeros(size(brain_mask.vol, 1), size(brain_mask.vol, 2), size(brain_mask.vol, 3), 'single');
  
  kernel = CBIG_gaussian3D(params.smooth_kernel_size, params.smooth_kernel_sigma);
  paddedKernel = CBIG_pad3DKernel(kernel, size(tmp_vol));
  fftK = fftn(paddedKernel);
  for d = 1:params.D
    word_indices              = params.word_indices{d};
    masked_vol(~word_indices) = 0;
    masked_vol(word_indices)  = 1;
    mask(brain_mask.vol(:) ~= 1) = 0;
    mask(brain_mask.vol(:) == 1) = masked_vol;

    smoothedMask = CBIG_FFTSmooth3DWithPredefinedKernel(mask, fftK);
    smoothedMask(mask ~= 1) = 0;
    for k = 1:params.K
      for a = 1:params.doc_by_num_authors(d)
        nonzero_vol(~word_indices) = 0;
        nonzero_vol(word_indices)  = squeeze(params.phi{d}(:, k, a));
        tmp_vol(brain_mask.vol(:) ~= 1) = 0;
        tmp_vol(brain_mask.vol(:) == 1) = nonzero_vol;
        
        tmp_vol = CBIG_FFTSmooth3DVolumeWithMaskPredfinedFFTKernel(tmp_vol, smoothedMask, fftK);
	    
        tmp = tmp_vol(brain_mask.vol(:) == 1);
        smoothed_phi{d}(:, k, a) = tmp(word_indices);
      end;
    end;
    denom = sum(sum(smoothed_phi{d},3),2);
    smoothed_phi{d} = bsxfun(@times, smoothed_phi{d}, 1./denom);
    smoothed_phi{d}(denom == 0,:,:) = 0;
  end;

  
  
function paddedKernel = CBIG_pad3DKernel(kernel, newSize)
  paddedKernel = zeros(newSize(1), newSize(2), newSize(3));
  kSize = size(kernel);
  paddedKernel(...
      uint64(newSize(1)/2-kSize(1)/2+1):uint64(newSize(2)/2+kSize(1)/2), ...
      uint64(newSize(2)/2-kSize(2)/2+1):uint64(newSize(2)/2+kSize(2)/2), ...
      uint64(newSize(3)/2-kSize(3)/2+1):uint64(newSize(3)/2+kSize(3)/2)) = kernel;

  
  
function h = CBIG_gaussian3D(dim, sigma)
  %   3D Gaussian lowpass filter
  %
  %   h = gausian3(dim,sigma) returns a rotationally
  %   symmetric 3D Gaussian lowpass filter with standard deviation
  %   sigma (in pixels). dim is a 1-by-3 vector specifying the number
  %   of rows, columns, pages in h. (dim can also be a scalar, in 
  %   which case h is dimxdimxdim.) If you do not specify the parameters,
  %   the default values of [3 3 3] for dim and 0.65 for
  %   sigma.


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
