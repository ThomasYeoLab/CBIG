function [ROI_PCs] = CBIG_preproc_aCompCor(func_file, ROI_file, nPCs, diff_flag)

% [ROI_PCs] = CBIG_preproc_aCompCor(func_file, ROI_file, nPCs)
% 
% Compute aCompCor regressors for a single volume. Detrend the regressors.
%
% Inputs:
%   - func_file:
%     Input BOLD filename (full path). E.g.
%     '<sub_dir>/Sub0001_Ses1/bold/002/Sub0001_Ses1_bld002_rest_skip4_stc_mc.nii.gz'
%
%   - ROI_file:
%     Filename of wm+ventricles mask (full path). E.g.
%     '<sub_dir>/Sub0001_Ses1/bold/mask/Sub0001_Ses1.func.wm.vent.nii.gz'
%
%   - nPCs:
%     Number of principle components extracted from wm+ventricles signals,
%     e.g. '5' or 5.
%
%   - diff_flag:
%     If diff_flag == 1, include first derivatives of components, otherwise
%     do not include derivatives.
%
% Outputs:
%   - ROI_PCs:
%     A T x nPCs matrix. The aCompCor regressors after detrending.
%
% Example:
%   [ROI_PCs] = CBIG_preproc_aCompCor('<sub_dir>/Sub0001_Ses1/bold/002/Sub0001_Ses1_bld002_rest_skip4_stc_mc.nii.gz', ...
%          '<sub_dir>/Sub0001_Ses1/bold/mask/Sub0001_Ses1.func.wm.vent.nii.gz', 5)
%
% Written by Jingwei Li.
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(ischar(nPCs))
    nPCs = str2num(nPCs);
end

if(ischar(diff_flag))
    diff_flag = str2num(diff_flag);
end

fmri = MRIread(func_file);
fmri_ROI = MRIread(ROI_file);

nframes = fmri.nframes;

if(sum(logical(fmri_ROI.vol(:))) == 0)
    ROI_PCs = nan(nframes, nPCs);
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick the voxels within the ROI mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vol = reshape(fmri.vol, size(fmri.vol, 1) * size(fmri.vol,2) * size(fmri.vol,3), nframes);
vol = vol(logical(fmri_ROI.vol), :)';               % TxN matrix
N = size(vol, 2);                                   % number of voxels in mask


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute mean of each voxel and the covariance across voxels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
covariance = zeros(nframes, nframes);
datamean = zeros(1, N);

% break up voxels to save memory
if(nframes <= 150)
    voxbinsize = 4000;
elseif(nframes > 150 && nframes <= 500)
    voxbinsize = 1000;
elseif(nframes > 500 && nframes <= 1000)
    voxbinsize = 500;
elseif(nframes > 1000 && nframes <= 3000)
    voxbinsize = 200;
else
    voxbinsize = 100;
end
for v = 1:voxbinsize:N
    idx = v:min(v + voxbinsize - 1, N);
    temp_vol = vol(:, idx);                         % T x voxbinsize matrix
    
    datamean(idx) = mean(temp_vol,1);
    temp_vol = bsxfun(@minus, temp_vol, datamean(idx));
    
    % compute covairance across voxels. If A = [A1 B1]; then A x A' = [A1
    % B1] x [A1 B1]' = A1 x A1' + B1 x B1'
    covariance = covariance + temp_vol * temp_vol';
end
clear temp_vol


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute principal components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[U, S, V] = svd(covariance);
dims = min(size(U, 2), nPCs);
PC = U(:, 1:dims) * diag(sqrt(diag(S(1:dims, 1:dims))));  
dPC = [zeros(1, dims); diff(PC)];

% output
if(diff_flag==1)
    ROI_PCs = [PC dPC];
else
    ROI_PCs = PC;
end
ROI_PCs = CBIG_glm_regress_matrix(ROI_PCs, [], 0);



end