function CBIG_MMLDA_brain_behavior_zscore_to_doc(Z_brain, Z_behavior, rid, ind_out, plus_num, nfold, out_dir, out_name)
% CBIG_MMLDA_brain_behavior_zscore_to_doc(Z_brain, Z_behavior, rid, ind_out, plus_num, nfold, out_dir, out_name)
%
% Convert brain behavior zscore to documents used to run MMLDA.
% 1. Plus a number to both brain and behavior zscore.
% 2. Threshold the zscore, we only consider non-negative values.
% 3. Times 10 to the zscore and discretize the zscore.
% 4. Output the subjects in <ind_out>
% 5. Split the output into <nfold> folds
%
% Input:
%   - Z_brain               : N x A marix, N is # of subjects, A is # of voxels within the mask.
%   - Z_behavior            : N x B matrix, N is # of subjects, B is # of scores.
%   - rid                   : cell array of id list 
%   - ind_out               : index of output cohort with respect to rid
%   - plus_num              : plus a number to z score of brain or behavior
%   - nfold                 : number of fold you want to split the output data (e.g., 1 or 10)
%   - out_dir               : output directory
%   - out_name              : output name of docs. For example, the output would be 
%                             If nfold == 1,
%                             <out_dir>/<out_name>_RID.txt
%                             <out_dir>/<out_name>_brain.dat
%                             <out_dir>/<out_name>_behavior.dat
%                             If nfold == 10,
%                             <out_dir>/<out_name>_RID_train<fold_num>.txt
%                             <out_dir>/<out_name>_RID_test<fold_num>.txt
%                             <out_dir>/<out_name>_brain_train<fold_num>.dat
%                             <out_dir>/<out_name>_brain_test<fold_num>.dat
%                             <out_dir>/<out_name>_behavior_train<fold_num>.dat
%                             <out_dir>/<out_name>_behavior_test<fold_num>.dat
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% plus a number to the z score, threshold and discretize
Z_brain = Z_brain+plus_num;
Z_brain(Z_brain < 0) = 0;
Z_brain = 10*Z_brain;
Z_brain = floor(Z_brain);

Z_behavior = Z_behavior+plus_num;
Z_behavior(Z_behavior < 0) = 0;
Z_behavior = 10*Z_behavior*size(Z_brain, 2)/size(Z_behavior, 2); 
Z_behavior = floor(Z_behavior);

% only output part of zscore
rid_out = rid(ind_out);
Z_brain_out = Z_brain(ind_out, :);
Z_behavior_out = Z_behavior(ind_out, :);

% divide z score into different folds and output to mmlda doc format
if nfold == 1
    csvwrite([out_dir '/' out_name '_RID.txt'], rid_out)
    CBIG_MMLDA_matrix2doc([out_dir '/' out_name '_brain.dat'], Z_brain_out)
    CBIG_MMLDA_matrix2doc([out_dir '/' out_name '_behavior.dat'], Z_behavior_out)
else 
    rng('default')
    rng(1);
    cv = cvpartition(size(Z_brain_out, 1), 'kfold', nfold);
    for i = 1:nfold
        csvwrite([out_dir '/' out_name '_RID_train' num2str(i) '.txt'], rid_out(cv.training(i)))
        csvwrite([out_dir '/' out_name '_RID_test' num2str(i) '.txt'], rid_out(cv.test(i)))
        CBIG_MMLDA_matrix2doc([out_dir '/' out_name '_brain_train' num2str(i) '.dat'], Z_brain_out(cv.training(i), :))
        CBIG_MMLDA_matrix2doc([out_dir '/' out_name '_behavior_train' num2str(i) '.dat'], Z_behavior_out(cv.training(i), :))
        CBIG_MMLDA_matrix2doc([out_dir '/' out_name '_brain_test' num2str(i) '.dat'], Z_brain_out(cv.test(i), :))
        CBIG_MMLDA_matrix2doc([out_dir '/' out_name '_behavior_test' num2str(i) '.dat'], Z_behavior_out(cv.test(i), :))
    end
end

