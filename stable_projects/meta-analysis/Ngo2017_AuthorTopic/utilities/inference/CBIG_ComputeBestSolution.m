function bestSeed = CBIG_ComputeBestSolution(root_dir, K, seeds, alpha, eta)
  % load mask
  mask_file = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2017_AuthorTopic', 'utilities', 'brain_mask', 'MNI_mask_conformed.2mm.0.1.nii.gz');
  mask = MRIread(mask_file);
  mask_index = find(mask.vol > 0);
  process_mat = zeros(length(mask_index), K, length(seeds));
  var_bound = zeros(length(seeds), 1);
  alphaEtaPrefix = ['alpha' num2str(alpha) '_eta' num2str(eta)];
  input_dir = fullfile(root_dir, ['K' num2str(K)], alphaEtaPrefix);

  % load samples
  disp('loading samples ...');
  tic
  count = 1;
  for seed = seeds
      disp(['seed: ' num2str(seed)]);
      load(fullfile(input_dir, ['params_K' num2str(K) '_SEED' num2str(seed) '.mat']));    
      if numel(mask_index) ~= size(params.beta, 2)
        error('Brain mask is different size from beta');
      end;
      process_mat(:, :, count) = params.beta';
      
      if(max(params.theta(:)) >= 1.1)
          var_bound(count) = -inf; 
      else
          var_bound(count) = params.var_bound(numel(params.var_bound));
      end
      count = count + 1;
  end
  toc

  mean(var_bound);
  [tmp, best_solution] = max(var_bound)


  disp('Computing correlation'); 
  tic
  process_mat_reshape = reshape(process_mat, [length(mask_index) K*length(seeds)]);
  z = CBIG_self_corr(process_mat_reshape);
  toc

  disp('Computing assignments to best solution');
  tic
  cost = nan(1, length(seeds));
  for i = best_solution
      for j = 1:length(seeds)
          if(i ~= j)
             row_start = (i - 1)*K + 1;
             row_end   = i*K;
             col_start = (j - 1)*K + 1;
             col_end   = j*K;
             x = z(row_start:row_end, col_start:col_end); 
              
             [tmp, cost(j)] = munkres(1 - x);  
          end
      end
  end
  toc

  % find solutions close to best solution
  disp('Finding solutions close to best');
  tic
  rel_cost = cost;
  rel_cost(best_solution) = 0;
  [closeness2best, solutions_close2best] = sort(rel_cost, 'ascend');

  % permute close solutions to match average solution
  matching_perm = zeros(length(seeds), K);
  for i = 1:length(seeds)
      [tmp, matching_perm(i, :)] = CBIG_corr_dist(process_mat(:, :, best_solution), process_mat(:, :, solutions_close2best(i)));
  end
  toc

  output_dir = fullfile(root_dir, 'best_solution', ['alpha' num2str(params.alpha) '_eta' num2str(params.eta)]);
  mkdir(output_dir);
  % save(fullfile(output_dir, ['BestSolution' num2str(K) '.mat']), 'best_solution', 'closeness2best', 'solutions_close2best', 'matching_perm');

  % save best solution
  seed = seeds(best_solution);
  load(fullfile(input_dir, ['params_K' num2str(K) '_SEED' num2str(seed) '.mat']));
  save(fullfile(output_dir, ['BestSolution_K' num2str(K, '%03d') '.mat']), 'params', 'closeness2best', 'solutions_close2best', 'matching_perm');
  bestSeed = seed;








      

% corr_dist
function [dist, matching] = CBIG_corr_dist(process1, process2)

  x = CBIG_my_corr(process1, process2);
  [matching, dist] = munkres(1 - x);
