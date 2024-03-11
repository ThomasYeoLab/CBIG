function [FC, FCD_CDF] = CBIG_pFIC_generate_GUSTO_subject_level_FC_FCD(TC_path, subject_id)
% [FC, FCD_CDF] = CBIG_pFIC_generate_GUSTO_subject_level_FC_FCD(TC_path, subject_id)
% This function generates a subject-specifc session-specific FC and FCD
% from the GUSTO dataset. Note that this function requires access to the subject's time
% series data, which are not open to public under the DUA. Thus this
% function is for reference and CBIG internal usage only.
%
% Input:
%   - TC_path: an absolute path to the directory containing subject time series
%   - subject_id: a char array of subject_id (for example, '21934' or '20620-012')
% Output:
%   - FC: a 68-by-68 subject-specific, session-specific FC
%   - FCD_CDF: a 10000-by-1 subject-specific FCD culumative distribution function (CDF)
%
% Example:
% [FC, FCD_CDF] = CBIG_pFIC_generate_GUSTO_subject_level_FC_FCD(['/isilon/CSC1/Yeolab/Data/GUSTO/pFIC/' ... 
%   'desikan_tc_noGSR/7.5'], '21215');
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% load parcellation time course
TC_path = dir([TC_path '/*' subject_id '_tc.mat']);
TC_filename = TC_path.name;
TC = load([TC_path '/' TC_filename]);
TC = TC.tc;
number_of_frame = size(TC, 2);
%% FC
FC = corr(TC');
%% FCD
TC = TC';
num_roi = size(FC, 1);
FC_size = num_roi*(num_roi-1)/2;
FCD_run = zeros(FC_size, number_of_frame - 22);
for j = 1:(number_of_frame - 22)
    TC_section = TC(j:j+22, :); 
    FC_section = CBIG_self_corr(TC_section); 
    FC_vec_section = FC_section(triu(true(size(FC_section, 1)), 1)); 
    FCD_run(:, j) = FC_vec_section;
end

FCD_run = corr(FCD_run); 
FCD_run_vec = FCD_run(triu(true(size(FCD_run, 1)), 1));
bin_count = histcounts(sort(FCD_run_vec), -1:0.0002:1);
FCD_CDF = cumsum(bin_count);

end
