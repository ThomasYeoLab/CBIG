function CBIG_AuthorTopicEM_RunGibbsSampler(seeds, K, ...
  burn_in, interval, num_samples, gibbs_input_path, output_dir, brain_mask, alpha, eta)
% CBIG_AuthorTopicEM_RunGibbsSampler(seeds, K, ...
%   burn_in, interval, num_samples, brain_mask, beta_init, theta_init)
%
% Estimate the author-topic model's parameters using Gibbs sampling.
%
% Input:
%  - seeds      : array of numbers denoting the indices of the
%                 initializations (seeds).
%  - K          : number of cognitive components.
%  - burn_in    : number burnin iteration of the Gibbs sampler
%  - gibbs_input_path: path to input file for the Gibbs sampler created
%                 by CBIG_AuthorTopicEM_ConvertActToGibbsFormat
%  - output_dir : folder containing GIBBS_outputs directory of the
%                 of Gibbs sampler's outputs
%  - brain_mask : struct containing the brain mask (read by MRIread)
%  - alpha & eta: hyperparameters of the author topic model's Dirichlet
%                 priors.
% Output:
%  Under <output_dir>/GIBBS_outputs/K<K>/alpha<alpha>_eta<eta>/seed<seed>/
%    - seed<seed>.burn<burn_in>.int<interval>.burnin.mat is the output from
%    the burn in
%    - seed<seed>.burn<burn_in>.sample<sample>.mat is the output from one of
%    the num_samples
%
% Example:
%   gibbs_input_path = '~/Work/gibbs_at_input.mat';
%   work_dir = '~/Work';
%   brain_mask2mm = MRIread(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
%     'meta-analysis', 'Yeo2015_AuthorTopicEM', 'utilities', ...
%     'mask', 'MNI_mask_conformed.2mm.0.1.nii.gz'));
%   CBIG_AuthorTopicEM_RunGibbsSampler(1:100, 2, 100, 100, 10, ...
%        gibbs_input_path, work_dir, brain_mask2mm, 100, 0.01);
%
%   Outputs of the Gibbs sampler are saved under
%   ~/Work/GIBBS_outputs/K2/alpha100_eta0.01/seed<seed> with seed between 1 and 100,
%   inclusive.
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
  if(nargin < 9)
     eta = 0.01;
     alpha = 100;
  end
  
  if(ischar(eta))
     eta = str2num(eta);
  end
  
  if(ischar(alpha))
     alpha = str2num(alpha);
  end
  
  if(ischar(seeds))
     seeds = str2num(seeds);
  end
  
  if(ischar(K))
     K = str2num(K);
  end
  
  if(ischar(burn_in))
     burn_in = str2num(burn_in);
  end
  
  if(ischar(interval))
     interval = str2num(interval);
  end
  
  if(ischar(num_samples))
     num_samples = str2num(num_samples);
  end
  
  % load mask
  mask_index = find(brain_mask.vol > 0);
  process_maps = brain_mask;
  process_maps.nframes = K;
  process_maps.vol = zeros([size(brain_mask.vol) K]);
  
  % load gibbs input
  load(gibbs_input_path);
  
  % Run LDA
  OUTPUT = 1;
  
  for seed = seeds
      seed_dir = fullfile(output_dir, 'GIBBS_outputs', ['K' num2str(K)], ...
        ['alpha' num2str(alpha) '_eta' num2str(eta)], ['seed' num2str(seed)]);
      system(['mkdir -p ' seed_dir]);
      
      % burn in
      tic
      disp('Running Gibbs sampler for burn in');
      disp(['  seed: ' num2str(seed)]);
      [WP, AT , Z , X ] = GibbsSamplerAT( WS , DS , paradigm_by_exp, K , burn_in , alpha , eta , seed , OUTPUT );
      
      phi = bsxfun(@times, full(WP) + eta , 1./(sum(full(WP), 1) + size(WP, 1) * eta));
      max_words = size(phi, 1);
      
      for i = 1:K
         brain_mask.vol(:) = 0;
         brain_mask.vol(mask_index(1:max_words)) = phi(:, i);
         process_maps.vol(:, :, :, i) = brain_mask.vol;
      end
      theta = bsxfun(@times, full(AT) + alpha, 1./(sum(full(AT), 2) + K * alpha));
      
      save(fullfile(seed_dir, ['seed' num2str(seed) ...
          '.burn' num2str(burn_in) '.int' num2str(interval) '.burnin.mat']), ...
          'process_maps', 'theta', 'alpha', 'eta');
  
      for sample = 1:num_samples
        [WP, AT , Z , X ] = GibbsSamplerAT( WS , DS , paradigm_by_exp, ...
          K , burn_in , alpha , eta , seed , OUTPUT , Z , X );
        
        phi = bsxfun(@times, full(WP) + eta , 1./(sum(full(WP), 1) + size(WP, 1) * eta));
        max_words = size(phi, 1);
        
        for i = 1:K
           brain_mask.vol(:) = 0;
           brain_mask.vol(mask_index(1:max_words)) = phi(:, i);
           process_maps.vol(:, :, :, i) = brain_mask.vol;
        end
        theta = bsxfun(@times, full(AT) + alpha, 1./(sum(full(AT), 2) + K * alpha));
        
      save(fullfile(seed_dir, ['seed' num2str(seed) '.burn' num2str(burn_in) ...
        '.int' num2str(interval) '.sample' num2str(sample) '.mat']), 'process_maps', 'theta', 'alpha', 'eta');
      end
      toc
  end
