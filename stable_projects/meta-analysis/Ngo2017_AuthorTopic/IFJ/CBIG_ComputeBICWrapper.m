function CBIG_ComputeBICWrapper(allKs)
  ALPHA = 100;
  ETA = 0.01;

  BIC_dir = 'BIC';
  smoothness_dir = fullfile(BIC_dir, 'ComponentSmoothness');
  mask_path = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2017_AuthorTopic', 'self-generated_thought', 'input_data', 'self-generated_thought_activation_mask.nii');

  component_dir = fullfile('outputs', 'best_solution', ['alpha' num2str(ALPHA) '_eta' num2str(ETA)]);
  solution_type = 'BestSolution';
  CBIG_EstComponentSmoothness(mask_path, allKs, component_dir, solution_type, BIC_dir);
  CBIG_ComputeBICAcrossK_byAFNISmoothness(allKs, component_dir, smoothness_dir, solution_type, BIC_dir);
