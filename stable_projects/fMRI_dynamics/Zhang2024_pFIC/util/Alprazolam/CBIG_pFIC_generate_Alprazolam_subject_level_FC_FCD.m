function [FC, FCD_CDF] = CBIG_pFIC_generate_Alprazolam_subject_level_FC_FCD(TC_path, subject_id, session, roi_list)

% [FC, FCD_CDF] = CBIG_pFIC_generate_Alprazolam_subject_level_FC_FCD(path, subject_id, session, roi_list)
% This function generates a subject-specifc session-specific FC and FCD from the Alprazolam
% dataset. Note that this function requires access to the subject's time
% series data, which are not open to public under the DUA. Thus this
% function is for reference and CBIG internal usage only.
%
% Input:
%   - TC_path: an absolute path to the directory containing subject time series
%   - subject_id: a 5-digit number (of 'double' type, not 'string')
%   - session: a binary value. 0 means drug session, 1 means placebo
%   session
%   - roi_list: a 72-by-1 binary vector, which entry reprensents if an ROI
%   is retained. 1 means the that ROI is retained, 0 means that the ROI is
%   excluded (due to ROI coverage <= 50% or medial wall)
% Output:
%   - FC: a 68-by-68 subject-specific, session-specific FC
%   - FCD_CDF: a 10000-by-1 subject-specific, session-specific FCD culumative distribution function (CDF)
%
% Example:
% [FC, FCD_CDF] = CBIG_pFIC_generate_Alprazolam_subject_level_FC_FCD(['isilon/CSC2/Yeolab/Data/Alpraz/' ... 
%   'TASK_fmriprep/desikan_from_vol/'], 13373, 0, 'replication/Alprazolam/input/drug/roi_list_0.5.txt')
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

TC = load([TC_path num2str(subject_id) '_' num2str(session) '_desikan_tc.mat']);
TC = TC.tc;
TC = TC';
fid = fopen(roi_list);
roi_list = textscan(fid, '%f');
fclose(fid);
roi_list = roi_list{1,1};
roi_list([1, 5, 37, 41]) = [];
num_roi = sum(roi_list, 1);
roi_list = repmat(roi_list, 1, 210);
TC = TC.*roi_list;
TC = TC(any(TC,2),:); % remove ROIs according to the ROI list, which represents ROIs with >50% coverage
%% FC
FC = corr(TC');
%% FCD
TC = TC';
FC_size = num_roi*(num_roi-1)/2;
FCD_run = zeros(FC_size, 210 - 19);
for j = 1:(210 - 19)
    TC_section = TC(j:j+19, :); 
    FC_section = CBIG_self_corr(TC_section); 
    FC_vec_section = FC_section(triu(true(size(FC_section, 1)), 1)); 
    FCD_run(:, j) = FC_vec_section;
end

FCD_run = corr(FCD_run); 
FCD_run_vec = FCD_run(triu(true(size(FCD_run, 1)), 1));
bin_count = histcounts(sort(FCD_run_vec), -1:0.0002:1);
FCD_CDF = cumsum(bin_count);
FCD_CDF = FCD_CDF';

end
