function homo_with_weight = CBIG_MSHBM_parameters_validation(project_dir,mesh,num_sess,num_clusters,subid,w,c)
% homo_with_weight = CBIG_MSHBM_parameters_validation(project_dir,mesh,num_sess,num_clusters,subid,w,c)
%
% This script will estimate the individual-level parcellation for the 
% validation subjects with num_clusters networks for a single subject 
% using num_sess sessions. The weight of group spatial prior should be 
% specified by w, the weight of MRF smoothness prior should be specified 
% by c.
% 
% After that, this script will apply the parcellation on the unseen
% validation fMRI files to compute the homogeneity <homo_with_weight>. 
% The validation fMRI files should be saved in 
% project_dir/data_list/validation_fMRI_list/lh_sub<?>.txt
% project_dir/data_list/validation_fMRI_list/lh_sub<?>.txt
% for data in 'fsaverage4/fsaverage5/fsaverage6/fsaverage'
% or
% project_dir/data_list/validation_fMRI_list/sub<?>.txt
% for data in 'fs_LR_32k'.
%
% The highter the <homo_with_weight>, the better the parcellation is. The
% user should run this sript with a grid search of w and c across all 
% subjects in validation set. And pick the pair of w and c with the highest
% <homo_with_weight>.
%
% To generate the individual-level parcellation, we assume the group priors
% are already estimated by running CBIG_MSHBM_estimate_group_priors.m. The 
% estimated group priors should be saved in 
% project_dir/priors/Params_Final.mat.
%
% We also assume the functional connectivity profiles of validation subjects 
% with num_sess sessions are already generated. The lists of functional 
% connectivity profiles should be saved in 
% project_dir/profile_list/validation_set/lh_sess<?>.txt
% project_dir/profile_list/validation_set/rh_sess<?>.txt
% for data in 'fsaverage4/fsaverage5/fsaverage6/fsaverage'
% or
% project_dir/profile_list/validation_set/sess<?>.txt
% for data in 'fs_LR_32k'.
%
% Input:
%   - project_dir:
%
%     The project directory.
%     1) project_dir/priors/Params_Final.mat.
%        contains the results of group priors estimated by running 
%        CBIG_MSHBM_estimate_group_priors.m on training set.
%        Params_Final.mat contain a struct variable Params with fields
%        corresponding to each group prior:
%        1) Params.epsil
%           Inter-subject functional connectivity variability. 
%        2) Params.mu
%           Group-level connectivity profiles for each network. 
%        3) Params.sigma
%           Intra-subject functional connectivity variability. 
%        4) Params.theta
%           Spatial prior which denotes the probability of each network 
%           occurring at each location. 
%                                                          
%     2) project_dir/profile_list/validation_set/lh_sess<?>.txt
%        project_dir/profile_list/validation_set/rh_sess<?>.txt
%        contain the functional connectivity profile lists of left
%        hemisphere and right hemisphere for each session. The functional
%        connectivity profiles are assumed to be pre-computed before run
%        the current script. These profile lists are assumed to be
%        pre-generated. For S test subjects and T sessions, there
%        should be T lh_sess<?>.txt and rh_sess<?>.txt lists for data in
%        'fsaverage4/fsaverage5/fsaverage6/fsaverage', or T sess<?>.txt
%        lists for data in 'fs_LR_32k'.
%        For example:
%        project_dir/profile_list/validation_set/lh_sess1.txt
%        project_dir/profile_list/validation_set/rh_sess1.txt
%        project_dir/profile_list/validation_set/lh_sess2.txt
%        project_dir/profile_list/validation_set/rh_sess2.txt
%        or
%        project_dir/profile_list/validation_set/sess1.txt
%        project_dir/profile_list/validation_set/sess2.txt
%        for 2 sessions. Each list should contain S rows, where each row 
%        is the full file path of the functional connectivity profile for 
%        each validation subject.
%
%     3) project_dir/data_list/validation_fMRI_list/lh_sub<?>.txt
%        project_dir/data_list/validation_fMRI_list/lh_sub<?>.txt
%        for data in 'fsaverage4/fsaverage5/fsaverage6/fsaverage'
%        or
%        project_dir/data_list/validation_fMRI_list/sub<?>.txt
%        for data in 'fs_LR_32k'
%        contain the validation fMRI files.The homogeneity will be computed
%        by applying parcellations on the validation fMRI files.
%        For example,if subject 1 has 2 sessions as validation fMRI files, 
%        session 1 has 1 fMRI file and session 2 has 2 fMRI files. Then
%        project_dir/data_list/validation_fMRI_list/lh_sub<?>.txt 
%        will contain:
%        lh_<fMRI_filename>_sub1_sess1_run1
%        lh_<fMRI_filename>_sub1_sess2_run1 lh_<fMRI_filename>_sub1_sess2_run2
%
%   - mesh: (string)
%     
%     The data surface space. 'fsaverage5/fsaverage6/fsaverage' or 'fs_LR_32k'. 
%
%   - num_sess: (string)
%
%     The number of sessions the user want to use to estimate the
%     individual-level parcellation. For example, '4'.
%
%   - num_clusters: (string)
%
%     The number of networks of the parcellations. For example, '17'.
%
%   - subid: (string)
%
%     The validation subject number. For example, '4' indicates the 4-th
%     subject in the project_dir/profile_list/validation_set/?h_sess<?>.txt or 
%     project_dir/profile_list/validation_set/sess<?>.txt.
%
%   - w: (string)
%   
%     The weight of group spatial prior Params.theta. For example, '100'.
%     A large w indicates strong weight of Params.theta. The estimated
%     individual-level parcellation will be very similar to the group-level
%     parcellation with very large w.
%
%   - c: (string)
%
%     The weight of MRF smoothness prior. For example, '50'. A large c
%     indicates more penalty for neighboring vertices being assigned to
%     different networks.
%
%
% Output:
%   
%   - homo_with_weight: 
%
%     The homogeneity value which is computed by applying the 
%     individual-level parcellation on the unseen fMRI files.The highter
%     the homogeneity, the better the parcellation quality. 
%     homo_with_weight is a vector with the same length as the validation 
%     fMRI list. Each fMRI file will has a homogeneity value. 
%
% Example:
%   homo_with_weight = CBIG_MSHBM_parameters_validation(project_dir,'fsaverage5','5','17','4','100','50')
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Generate parcellations for the validation set
[lh_labels, rh_labels] = CBIG_MSHBM_generate_individual_parcellation( ...
                         project_dir, mesh, num_sess, num_clusters, subid, w, c, 'validation_set');

%% Compute homogeneity for the validation set

if(~isempty(strfind(mesh,'fsaverage')))
    lh_fMRI_files = fullfile(project_dir,'data_list','validation_fMRI_list',['lh_sub' subid '.txt']);
    rh_fMRI_files = fullfile(project_dir,'data_list','validation_fMRI_list',['rh_sub' subid '.txt']);
    homo_with_weight = CBIG_ComputeParcellationHomogeneity_FS(lh_labels,rh_labels,mesh,lh_fMRI_files,rh_fMRI_files);
    
elseif(~isempty(strfind(mesh,'fs_LR_32k')))
    fMRI_files = fullfile(project_dir,'data_list','validation_fMRI_list',['rh_sub' subid '.txt']);
    homo_with_weight = CBIG_ComputeParcellationHomogeneity_fslr(lh_labels,rh_labels,fMRI_files);
end

%% Save results
out_dir=fullfile(project_dir,'homogeneity', 'validation_set');

if(~exist(out_dir))
    mkdir(out_dir);
end
save(fullfile(out_dir, ...
    ['Ind_homogeneity_MSHBM_sub',num2str(subid),'_w',num2str(w),'_MRF',num2str(c),'.mat']), ...
    'homo_with_weight');
