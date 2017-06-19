function CBIG_EstComponentSmoothness(mask_path, allKs, component_dir, solution_type, output_dir)
  VOXEL_DIM = 2; % mm
  mask_brain = MRIread(mask_path);
  spm_mask_brain = spm_vol(mask_path);
  search_vol = sum(mask_brain.vol(:)) * VOXEL_DIM^3;
  
  output_dir = fullfile(output_dir, 'ComponentSmoothness');
  system(['mkdir -p ' output_dir]);
  
  full_brain_mask = MRIread('~ngohgia/templates/MNI_mask_conformed.2mm.0.1.nii.gz');
  % Smoothness estimation
  for K = allKs
    solution = load(fullfile(component_dir, [solution_type '_K' num2str(K, '%03d')]));
    if exist(fullfile(output_dir, [solution_type '_AFNI_K' num2str(K) '_smoothness.txt']), 'file')
      system(['rm ' fullfile(output_dir, [solution_type '_AFNI_K' num2str(K) '_smoothness.txt'])]);
    end
    
    for C = 1:K
      [component_path, znorm_component_path] = CBIG_SaveATComponent(solution.params, output_dir, full_brain_mask, solution_type);
      % Get AFNI's FWHM estimation
      system(['3dFWHMx -input ' znorm_component_path ' -mask ' mask_path ' >> ' fullfile(output_dir, [solution_type '_AFNI_K' num2str(K) '_smoothness.txt'])]);
    end
  end
  
  % perform smoothness estimation
  for K = allKs  
    afni_smoothness_file = fullfile(output_dir, [solution_type '_AFNI_K' num2str(K) '_smoothness.txt']);
    fileID = fopen(afni_smoothness_file, 'r');
    
    afni_fwhm = zeros(K, 3);
    spm_fwhm = zeros(K, 3);
    afni_resels_count = zeros(K, 1);
    afni_spm_resels_count = zeros(K, 1);
    spm_resels_count = zeros(K, 1);
    
    for i = 1:K
      % using AFNI's FWHM estimation and self-computed resels count
      line = fgetl(fileID);
      curr_fwhm = strread(line, '%f');
      afni_fwhm(i, :) = curr_fwhm';
      fwhm_vol = curr_fwhm(1) * curr_fwhm(2) * curr_fwhm(3);
      afni_resels_count(i) = search_vol / fwhm_vol;
    end
    sum_afni_resels_count = sum(afni_resels_count);
    save(fullfile(output_dir, ['K' num2str(K) '_smoothness.mat']), 'afni_fwhm', 'sum_afni_resels_count');
    fclose(fileID);
  end
