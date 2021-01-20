function flags = CBIG_MMLDA_check_example_results(out_dir)
%
% This function compares the example results saved in <out_dir> with the
% reference results.
%
% Input:
%     - out_dir:
%       The output directory where example results are saved.
%
% Output:
%     - flags:
%       flags is logical 5x1 vector that indicates whether the example results
%       saved in <output_dir> match the reference results. Vector is in the 
%       order of beta1, beta2, gamma, likelihood, inferred gamma. If example 
%       results match reference results within tolerance, flag = 1, otherwise 
%       flag = 0.
%
% Written by Leon OOI and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    % define output folder
    CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
    % check the results
    ref_dir = fullfile(CBIG_CODE_DIR,'stable_projects/disorder_subtypes', ...
    'Sun2019_ADJointFactors/examples/correct_output');
    flags = logical(ones(1,5));
    % compare beta1, beta2 and gamma estimates
    params = {'beta1' 'beta2' 'gamma'};
    curr_out_dir = fullfile(out_dir, 'estimation/ADNI2_bl_AD_meanCNstdALL_plus1_2sub/k2/r3');
    curr_ref_dir = fullfile(ref_dir, 'estimation/ADNI2_bl_AD_meanCNstdALL_plus1_2sub/k2/r3');
    for n = 1:length(params)
        est = load(fullfile(curr_out_dir, strcat('final.', params{n})));
        ref_est = load(fullfile(curr_ref_dir, strcat('final.', params{n})));
        rel_diff = max(max(abs((est - ref_est) ./ ref_est)));
        if rel_diff > 1e-4; 
            sprintf('maximum %s rel difference: %f', params{n}, rel_diff)
            flags(n) = 0;
        end
    end

    % compare likelihood estimate
    likelihood = load(fullfile(curr_out_dir, 'likelihood.dat'));
    likelihood = likelihood(:, 1);
    ref_likelihood = load(fullfile(curr_ref_dir, 'likelihood.dat'));
    ref_likelihood = ref_likelihood(:, 1);
    rel_diff_likelihood = max(max(abs((likelihood - ref_likelihood) ./ ref_likelihood)));
    if rel_diff > 1e-4; 
            sprintf('maximum likelihood rel difference: %f', rel_diff_likelihood)
            flags(4) = 0;
    end

    % compare inferred gamma file
    curr_out_dir = fullfile(out_dir, 'inference/ADNI2_bl_AD_meanCNstdALL_plus1_2sub');
    curr_ref_dir = fullfile(ref_dir, 'inference/ADNI2_bl_AD_meanCNstdALL_plus1_2sub');
    gamma = load(fullfile(curr_out_dir, 'k2_inf_ADNI2_bl_AD_meanCNstdALL_plus1_1sub-gamma.dat'));
    ref_gamma = load(fullfile(curr_ref_dir, 'k2_inf_ADNI2_bl_AD_meanCNstdALL_plus1_1sub-gamma.dat'));
    rel_diff_gamma = max(max(abs((gamma - ref_gamma) ./ ref_gamma)));
    if rel_diff > 1e-4; 
            sprintf('maximum inferred gamma rel difference: %f', rel_diff_gamma)
            flags(5) = 0;
    end
end
