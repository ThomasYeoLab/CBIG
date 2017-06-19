function [component_path, znorm_component_path] = SaveComponent(params, output_dir, brain_mask, prefix)
  if nargin == 2
    prefix = '';
  end
  
  if nargin == 1
    brain_mask = MRIread('~ngohgia/templates/MNI_mask_conformed.2mm.0.1.nii.gz');
  end
  
  K = size(params.beta, 1);
  
  for i = 1:K
    new_brain = brain_mask;
    new_brain.vol = zeros(size(new_brain.vol));
    vol = params.beta(i, :);
    new_brain.vol(brain_mask.vol ~= 0) = vol;
    component_path = fullfile(output_dir, [prefix '_K' num2str(K) '_C' num2str(i) '.nii']);
    MRIwrite(new_brain, component_path);
    
    znorm_vol = zscore(vol);
    new_brain.vol = zeros(size(new_brain.vol));
    new_brain.vol(brain_mask.vol ~= 0) = znorm_vol;
    znorm_component_path = fullfile(output_dir, [prefix '_K' num2str(K) '_C' num2str(i) '_znorm.nii']);
    MRIwrite(new_brain, znorm_component_path);
  end