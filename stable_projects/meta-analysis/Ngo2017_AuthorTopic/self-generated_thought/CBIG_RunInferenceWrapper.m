function CBIG_RunInferenceWrapper(allKs, seeds)
  ALPHA     = 100;
  ETA       = 0.01;
  SMOOTH    =1;
  
  BASE_DIR  = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2017_AuthorTopic', 'self-generated_thought');
  DATA_PATH = fullfile(BASE_DIR, 'input_data', 'self-generated_thought.mat');

  for K = allKs
    for seed = seeds
      CBIG_RunAT_CVB(seed, K, ALPHA, ETA, SMOOTH, BASE_DIR, DATA_PATH);
    end
  end
