function score_pca = CBIG_pFIC_PCA_GUSTO_behavior_score(TC_path, subject_list)

% score_pca = CBIG_pFIC_PCA_GUSTO_behavior_score(TC_path, subject_list)
% This function computes the first PCA component of the 5 behavior scores from GUSTO dataset
% that are used in the manuscript. The 5 behavior scores include:
% (1) Cambridge Neuropsychological Test Automated Battery (CANTAB) Spatial Working Memory (SWM) test 
%   (completed at age 6): sum of total errors for 4 and 6 boxes trails.
% (2) Delayed Matching to Sample (DMS) test (completed at age 6): percentage of the total number of trials 
%   upon which a correct selection was made on the participant's first response.
% (3) Behavior Rating Inventory of Executive Function (BRIEF; completed at age 7): Cognition Regulation Index T-score.
% (4) Wechsler Abbreviated Scale of Intelligence (WASI) test (completed at age 7): sum of Block Design and Matrix 
%   Reasoning T-scores.
% (5) CANTAB SWM test (completed at age 8.5): sum of total errors for 4 and 8 boxes trails.
%
% Input:
%   - TC_path: an absolute path to the directory containing behavioral scores
%   - subject_list: an absolute path to the list of subjects
% Output:
%   - score_pca: 1st PC of the 5 behavior scores mentioned above
%
% Example:
% score_pca = CBIG_pFIC_PCA_GUSTO_behavior_score(['/isilon/CSC1/Yeolab/Data/GUSTO/pFIC/behavior/' ... 
%   'clean_scores/'], '/home/shaoshi.z/storage/MFM/GUSTO/behavior/subject_list/subject_list_sorted_age.txt')
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Define the list of behvaior scores considered
behavior_score_list = {[TC_path '/Y06_SWM.csv'],...
    [TC_path '/Y06_DMS.csv'], ...
    [TC_path '/Y07_BRIEF.csv'], ...
    [TC_path '/Y07_WASI.csv'], ...
    [TC_path '/Y08.5_SWM.csv']}; 
score = [];  

%% loop through the list of behavior scores and extract the relevant scores
for i = 1:length(behavior_score_list)
    behavior = readtable(behavior_score_list{i});
    if strcmp(behavior_score_list{i}, [TC_path '/Y06_DMS.csv'])
        ind = extract_behavior(subject_list, behavior.ID);
        subscore = behavior.DMS_Percent_correct_yr6(ind);
    elseif strcmp(behavior_score_list{i}, [TC_path '/Y06_SWM.csv'])
        ind = extract_behavior(subject_list, behavior.ID);
        subscore = behavior.SWM_Total_errors_4_to_6_boxes_yr6(ind);
    elseif strcmp(behavior_score_list{i}, [TC_path '/Y07_WASI.csv'])
        ind = extract_behavior(subject_list, behavior.ID);
        subscore = [behavior.PerceptualReasoningSumOfT_scores(ind)];
    elseif strcmp(behavior_score_list{i}, [TC_path '/Y07_BRIEF.csv'])
        ind = extract_behavior(subject_list, behavior.ID);
        subscore = behavior.brief2pf_cri_t(ind);
    elseif strcmp(behavior_score_list{i}, [TC_path '/Y08.5_SWM.csv'])
        ind = extract_behavior(subject_list, behavior.ID);
        subscore = behavior.SWMTotalErrors_4To8Boxes_(ind); 
    end
    score = [score subscore];
end

%% normalize the scores and extract the first principal component
score_z = normalize(score);
[~, score_pca] = pca(score_z);
score_pca = score_pca(:, 1);

end

%% helper function
% This function loads a list of subject IDs and find a match from the
% behavior score list (the ID column specifically)
function ind = extract_behavior(subject_list, behavior_sub_list)
fid = fopen(subject_list);
subject_list = textscan(fid, '%s');
subject_list = subject_list{1};
fclose(fid);

ind = [];
for i = 1:length(subject_list)
    ind(i) = find(contains(behavior_sub_list, subject_list{i}));
end
end