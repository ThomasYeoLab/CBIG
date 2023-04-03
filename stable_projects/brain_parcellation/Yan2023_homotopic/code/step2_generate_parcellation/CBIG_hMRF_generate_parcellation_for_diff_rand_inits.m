function [params, results] = CBIG_hMRF_generate_parcellation_for_diff_rand_inits(params)
% [params, results] = CBIG_hMRF_generate_parcellation_for_diff_rand_inits(params)
%
% This function runs the clustering algorithm based on the input parameter structure.
% 
% Input
%   - params: (struct) 
%     The structure that consists of all the customized parameters that the user would like to pass in.
%     See CBIG_hMRF_set_params.m for more details about what parameters the user could set.
%
% Example
%   - [params, results] = CBIG_hMRF_generate_parcellation_for_diff_rand_inits(params)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
disp('Loading PMM matrix...');
lhrh_vol = load(params.lhrh_avg_file);

disp('##########################');
disp('# Random initializations #');
disp('##########################');

results_dir = fullfile(params.output_folder, 'results');
if(~exist(results_dir, 'dir'))
    mkdir(results_dir);
end

for i = params.start_seed:(params.start_seed + params.num_rand_inits - 1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set up seed-specific parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.seed = i;
    
    % the debug_out_folder can be used to save any temporary files specific to the current seed
    debug_out_folder = fullfile(params.output_folder, ['seed_' num2str(params.seed)]);
    if(~exist(debug_out_folder, 'dir'))
        mkdir(debug_out_folder);
    end
    params.debug_out_folder = debug_out_folder;

    % set up convergence file path
    params.convergence_log_path = fullfile(params.convergence_log_folder,...
     ['convergence_log_' num2str(params.seed) '.txt']);

    fprintf('This is using seed %d \n.', i);
    if(exist([params.output_folder, params.output_name, '_seed_', num2str(params.seed), '_Energy.mat'],'file'))
        continue;
    else
        disp('###########################################');
        disp('# Perform clustering for the current seed #');
        disp('###########################################');
        
        results = CBIG_hMRF_optimize_cost_function(lhrh_vol.final_PMM, lhrh_vol.dim, params);
        save(fullfile(results_dir, [params.output_name, '_seed_', num2str(params.seed), '.mat']),...
            'results', 'params');
        fprintf('Result for seed %d has been saved.\n', params.seed);
    end
end
fprintf('SUCCESS!\n');

end
