function CBIG_AuthorTopic_EstimateComponentSmoothness(allKs, maskPath, ...
   bestSolutionDir, outputDir)
% CBIG_AuthorTopic_EstimateComponentSmoothness(allKs, maskPath, ..
%   bestSolutionDir, outputDir)
%
% Estimate the smoothness of the estimates for the probability of a
% component activating a voxel (Pr(voxel | component)), at different
% numbers of cognitive components.
% The smoothness is estimated as the FWHM of a 3D Gaussian smoothing kernel
% if the kernel were used to smooth the model parameter estimates.
%
% Input:
%  - allKs          : vector containing the numbers of components.
%  - maskPath       : path to the mask defining the voxels that are used
%                     for BIC computation. This should be the mask of
%                     activation across all experiments included in the
%                     meta-analysis, produced by
%                     CBIG_AuthorTopic_GenerateCVBDataFromText.m
%  - bestSolutionDir: path to directory containing the best model
%                     parameters' estimates for each given number of
%                     cognitive components.
%  - outputDir      : directory containing all intermediate files.
% Output:
%
% Example:
%   CBIG_AuthorTopic_EstimateComponentSmoothness([1:5], ...
%     '/Work/AT/expMask.nii.gz', '/Work/AT/outputs/bestSolution', ...
%     '/Work/AT/BIC')
%   Estimate the smoothness of the author-topic model's probability of a
%   component activating a voxel (Pr(voxel | component)) for 1 to 5
%   cognitive components. The estimates are from
%   '/Work/AT/outputs/bestSolution/BestSolution_K001.mat' to
%   '/Work/AT/outputs/bestSolution/BestSolution_K005.mat'. The mask comes
%   from '/Work/AT/expMask.nii.gz'. All output files are saved under
%   '/Work/AT/BIC'.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  % get voxel dimension from one of the solutions
  firstK = allKs(1);
  firstSolution = load(fullfile(bestSolutionDir, ['BestSolution_K' num2str(firstK, '%03d')]));
  fullBrainMask = MRIread(firstSolution.params.maskPath);

  voxelDim = fullBrainMask.volres(1);

  maskBrain = MRIread(maskPath);
  searchVol = sum(maskBrain.vol(:)) * voxelDim^3;

  % Estimate FWHM of the smoothness in component's brain image with AFNI
  for K = allKs
    solution = load(fullfile(bestSolutionDir, ['BestSolution_K' num2str(K, '%03d')]));
    afniOutputFile = fullfile(outputDir, ['BestSolution_AFNI_K' num2str(K) '_smoothness.txt']);
    if exist(afniOutputFile, 'file')
      system(['rm ' afniOutputFile]);
    end

    for C = 1:K
      znormBrain = fullBrainMask;
      znormBrain.vol = zeros(size(znormBrain.vol));

      znormBeta = zscore(solution.params.beta(C, :));
      znormBrain.vol(fullBrainMask.vol ~= 0) = znormBeta;
      znormBrainImagePath = fullfile(outputDir, ['BestSolution_K' num2str(K) '_C' num2str(C) '_znorm.nii']);
      MRIwrite(znormBrain, znormBrainImagePath);

      % Call AFNI's 3dFWHMx
      system(['3dFWHMx -input ' znormBrainImagePath ' -mask ' maskPath ' >> ' afniOutputFile]);
    end
  end

  % read AFNI's 3dFWHMx output
  for K = allKs
    afniOutputFile = fullfile(outputDir, ['BestSolution_AFNI_K' num2str(K) '_smoothness.txt']);
    fileID = fopen(afniOutputFile, 'r');

    afniFWHM = zeros(K, 3);
    reselsCount = zeros(K, 1);

    for C = 1:K
      line = fgetl(fileID);
      currFwhm = strread(line, '%f');
      afniFWHM(C, :) = currFwhm';
    end
    minAfniFWHM = min(afniFWHM, [], 1);
    fwhmVol = minAfniFWHM(1) * minAfniFWHM(2) * minAfniFWHM(3);
    totalReselsCount = searchVol / fwhmVol * K;
    save(fullfile(outputDir, ['K' num2str(K) '_smoothness.mat']), 'afniFWHM', 'totalReselsCount');
    fclose(fileID);
  end
