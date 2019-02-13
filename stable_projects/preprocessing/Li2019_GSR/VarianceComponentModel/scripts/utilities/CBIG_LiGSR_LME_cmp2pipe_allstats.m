function CBIG_LiGSR_LME_cmp2pipe_allstats( trait_list, pipeline1_dir, pipeline2_dir, num_families, num_samples, d, outdir )

% CBIG_LiGSR_LME_cmp2pipe_allstats( trait_list, pipeline1_dir, pipeline2_dir, num_families, num_samples, d, outdir )
% 
% This function calculates all statistics mentioned in Li et al., 2019 when
% comparing the two preprocessing pipelines using variance component model,
% including the percentage improvement of explained variance, the number of
% traits with the entire inter-quartile range (IQR) of the explained
% variance difference above zero (or below zero), the number of traits
% with the median of the explained variance above zero (or below zero).
% 
% Inputs:
%   - trait_list
%     A string. The full path of the list of all traits involved in the
%     variance component model estimates.
% 
%   - pipeline1_dir
%     A string. The variance component model output directory for the first
%     preprocessing pipeline. It is assumed there are a bunch of subfolders
%     named as ['del' num2str(d) '_set' num2str(l)] for l-th delete-d
%     jackknife sample. 
%     This function will read the following .mat files:
%     fullfile(pipeline1_dir, ['del' num2str(d) '_set' num2str(l)], ...
%         ['m2' quant_stem '_' traits{i} '.mat'])
%     where "d" is the number of subject IDs removed for each jackknife
%     sample (see below); "l" ranges from 1 to num_samples (see
%     below); quant_stem = '' when the trait name is Sex or Gender,
%     otherwise quant_stem = '_QuantileNorm'; i ranges from 1 to the total
%     number of traits.
% 
%   - pipeline2_dir
%     A string. The variance component model output directory for the second
%     preprocessing pipeline. It is assumed there are a bunch of subfolders
%     named as ['del' num2str(d) '_set' num2str(l)] for l-th delete-d
%     jackknife sample.
%     This function will read the following .mat files:
%     fullfile(pipeline2_dir, ['del' num2str(d) '_set' num2str(l)], ...
%         ['m2' quant_stem '_' traits{i} '.mat'])
%     where "d" is the number of subject IDs removed for each jackknife
%     sample (see below); "l" ranges from 1 to num_samples (see
%     below); quant_stem = '' when the trait name is Sex or Gender,
%     otherwise quant_stem = '_QuantileNorm'; i ranges from 1 to the total
%     number of traits.
% 
%   - num_families
%     A scalar, the number of families in the full subject list. If all
%     subjects are unrelated, then the number of families equals to the
%     number of subjects.
% 
%   - num_samples
%     A scalar, the number of jackknife samples.
%  
%   - d
%     A scalar, the number of subjects that were removed for each jackknife
%     sample.
% 
%   - outdir
%     A string. The full path of the output directory. A file
%     "allstats_cmp2pipelines.mat" will be created in the output directory.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

[traits, num_trait] = CBIG_text2cell(trait_list);

for i = 1:num_trait
    for sample = 1:num_samples
        if (strcmpi(traits{i}, 'Sex') || strcmpi(traits{i}, 'Gender'))
            quant_stem = '';
        else
            quant_stem = '_QuantileNorm';
        end
        pipe1_name = fullfile(pipeline1_dir, ['del' num2str(d) '_set' num2str(sample)], ...
            ['m2' quant_stem '_' traits{i} '.mat']);
        pipe2_name = fullfile(pipeline2_dir, ['del' num2str(d) '_set' num2str(sample)], ...
            ['m2' quant_stem '_' traits{i} '.mat']);
        
        load(pipe1_name)
        m_del_d_mat1(sample, i) = morpho.m2;    % num_samples x num_trait
        clear morpho
        load(pipe2_name)
        m_del_d_mat2(sample, i) = morpho.m2;    % num_samples x num_trait
        clear morpho
    end
end

[perc_improv, m_jack, v_jack] = CBIG_LiGSR_del_d_jack_cmp2pipelines( ...
    m_del_d_mat1, m_del_d_mat2, num_families, d );
[IQR_pos, IQR_neg, med_pos, med_neg, mean_m2] = CBIG_LiGSR_PosNeg_jackIQR_cmp2pipe( ...
    m_del_d_mat1, m_del_d_mat2 );

if(~exist(outdir, 'dir'))
    mkdir(outdir)
end
save(fullfile(outdir, 'allstats_cmp2pipelines.mat'), 'perc_improv', 'm_jack', 'v_jack', ...
    'IQR_pos', 'IQR_neg', 'med_pos', 'med_neg', 'mean_m2')

end

