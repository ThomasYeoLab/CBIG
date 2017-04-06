function CBIG_brain2doc(vol4D, concatOrder, mask, refListOrParams, nuisanceVars, outDir, varargin)

% CBIG_brain2doc(vol4D, concatOrder, mask, refListOrParams, nuisanceVars, outDir, varargin)
%
% This function converts registered, modulated GM images into documents
% that can be analyzed by LDA. It will first take log10 of the voxel
% values, regress out nuisance variables for all images w.r.t. the
% "reference group", z-normalize each voxel of all images w.r.t. the same
% reference group, and finally output the document files, which will later
% serve as the input of LDA.
%
% Input:
%     - vol4D:
%       MRI 4D volume
%     - concatOrder:
%       Text file specifying the concatenation order of vol4D
%       This file is essentially a reordered version of nameList, input of
%       ../../step1_VBM/extractBrains/lib/CBIG_runBET.sh
%     - mask: 
%       Binary GM mask (MRI volume)
%     - refListOrParams: 
%       Text file or .mat file
%       If the former, it is a subset of concatOrder, specifying the
%       "reference group", w.r.t. which regression and z-normalization are
%       performed. If the latter, the function will load the regression and
%       z-normalization parameters from it (instead of computing with the
%       reference group)
%     - nuisanceVars: 
%       .csv file with each line specifying the nuisance variables for the
%       corresponding image in concatOrder
%       The effects of these nuisance variables will be regressed out by
%       GLM
%     - outDir:
%       Output directory
%     - (optional) filenameList1, filenameList2, ...
%       Text file specifying images that you want to cluster into one file
%       Therefore, it is a subset of concatOrder. For example, if you want
%       two files with the first containing word counts of images 1, 2 and
%       3 and the second containing images 1, 2 and 4, just input two text
%       files specifying so. If not provided, all images will be written
%       into a single file
%
% Output:
%     - .dat file(s) summarizing word counts of the images 
%       The number output files generated is determined by how many
%       filenameList arguments you provided. If at least one filenameList is
%	provided, the image order follows filenameList(s). Otherwise, it
%	follows concatOrder
%   
% Example:
% See ../replicatePNAS/CBIG_LDA_wrapper.m
%
% Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Load data
% Load 4D MRI
mri = MRIread(vol4D);
vols = mri.vol;
% Load image order
fID = fopen(concatOrder);
imgOrder = textscan(fID, '%s\n');
fclose(fID);
imgOrder = imgOrder{1};
% Load mask
mri = MRIread(mask);
mask_3d = mri.vol;
% Load reference list or GLM parameters
if strcmp(refListOrParams(end-3:end), '.mat') % GLM parameters given
    paramGiven = true;
    load(refListOrParams);
else % reference list given
    paramGiven = false;
    fID = fopen(refListOrParams);
    refs = textscan(fID, '%s\n');
    fclose(fID);
    refs = refs{1};
end
% Load nuisance regressors
regressors = csvread(nuisanceVars);
warning('Rows of concatOrder and nuisanceVars must correspond. If not, this function will not work properly.');


%% Apply mask and take log10
noImg = size(vols, 4);
noVox = sum(mask_3d(:));
volMat = zeros(noImg, noVox); % each row holds an image's all voxels
for idx = 1:noImg
    vol_3d = vols(:, :, :, idx);
    % voxel coordinates (i, j, k) in FreeView corresponds to
    % matrix coordinates (j+1, i+1, k+1) in MATLAB
    vol_1d = reshape(vol_3d, [1 numel(vol_3d)]);
    mask_1d = reshape(mask_3d, [1 numel(mask_3d)]);
    vol_1d(mask_1d==0) = []; % remove voxels not in the GM mask
    vol_1d(vol_1d==0) = 1; % set voxels in non-GM area to 1, so that later log(1) = 0, i.e., no atrophy
    volMat(idx, :) = log(vol_1d)/log(10); % take log10()
end
disp('Masking and Taking log10 -- Finished');


%% Regress out nuisance variables
% Estimating GLM parameters with only images in refList
if ~paramGiven
    refInd = ismember(imgOrder, refs);
    refY = volMat(refInd, :);
    refYMean = mean(refY, 1);
    refX = regressors(refInd, :);
    % GLM fitting column by colume (voxel by voxel)
    refGLMBetas = cell(noVox, 1);
    for idx = 1:noVox
        [betas, ~, ~] = glmfit(refX, refY(:, idx));
        refGLMBetas{idx} = betas;
    end
    disp('GLM Parameter Estimation -- Finished');
end
% Regress out nuisance variables for ALL images by computing
% y - X*betas_ref + mean(y_ref) for each voxel
X = [ones(size(regressors, 1), 1) regressors];
volMat_reg = zeros(noImg, noVox);
for idx = 1:noVox % for each voxel
    y = volMat(:, idx);
    betas = refGLMBetas{idx};
    volMat_reg(:, idx) = y-X*betas+refYMean(idx);
end
disp('Regression -- Finished');


%%  Z-score and discretize
% Compute post-regression reference group mean and std
mkdir(outDir);
if ~paramGiven
    refYMean_reg = mean(volMat_reg(refInd, :), 1);
    refYStd_reg = std(volMat_reg(refInd, :), 0, 1);
    save([outDir 'refParams.mat'], 'refYMean', 'refGLMBetas', 'refYMean_reg', 'refYStd_reg');
end
% Z-normalization for each voxel using reference group's mean and std
Z = zeros(size(volMat_reg));
for idx = 1:noVox
    m = refYMean_reg(idx);
    s = refYStd_reg(idx);
    Z(:, idx) = (volMat_reg(:, idx)-m)/s;
end
Z(Z>0) = 0; % higher than mean is set to 0 (considered as no atrophy)
tc = floor(-10*Z); % tc for "term (dictionary word) counts"

%% Output documents
disp('Writing to files...');
if numel(varargin) == 0
    % Output a single document containing all images
    fID = fopen([outDir 'docs_all.dat'], 'w');
    fID = fopen([outDir 'docs_all.dat'], 'a');
    for idx1 = 1:noImg
        tc_oneImg = tc(idx1, :);
        noTerms = sum(tc_oneImg~=0);
        fprintf(fID, '%i ', noTerms);
        for idx2 = 1:numel(tc_oneImg)
            if tc_oneImg(idx2) ~= 0
                fprintf(fID, '%i:%i ', idx2-1, tc_oneImg(idx2));
            end
        end
        fprintf(fID, '\n');
    end
    fclose(fID);
else
    % Output one or multiple documents as requested
    for idx1 = 1:numel(varargin)
        % Read list
        filenameList = varargin{idx1};
        fID = fopen(filenameList);
        filenames = textscan(fID, '%s\n');
        fclose(fID);
        filenames = filenames{1};
        % Fetch term counts for each filename
        fID = fopen([outDir 'docs_' inputname(6+idx1) '.dat'], 'w');
        fID = fopen([outDir 'docs_' inputname(6+idx1) '.dat'], 'a');
        for idx2 = 1:numel(filenames)
            logIdx = strcmp(imgOrder, filenames{idx2});
            tc_oneImg = tc(logIdx, :);
            noTerms = sum(tc_oneImg~=0);
            fprintf(fID, '%i ', noTerms);
            for idx3 = 1:numel(tc_oneImg)
                if tc_oneImg(idx3) ~= 0
                    fprintf(fID, '%i:%i ', idx3-1, tc_oneImg(idx3));
                end
            end
            fprintf(fID, '\n');
        end
        fclose(fID);
    end
end
