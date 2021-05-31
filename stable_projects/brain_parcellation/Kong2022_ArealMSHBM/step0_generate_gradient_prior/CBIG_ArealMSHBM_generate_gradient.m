function CBIG_ArealMSHBM_generate_gradient(mesh,out_dir,sub,num_sess)

% CBIG_ArealMSHBM_generate_gradient(mesh,out_dir,sub,num_sess)
% This function will generate diffusion embedding matrices of
% gradients for a specific subject with <num_sess> session.  
%
% Input:
%   - mesh:
%
%     The surface space of fMRI data, e.g. 'fsaverage5' or 'fs_LR_32k'.
%     The data is allowed to be in either fsaverage space (e.g. 
%     fsaverage4/5/6, fsaverage) or fs_LR_32k. Note that fs_LR_164k is not
%     available.
%
%   - out_dir:
%
%     1) fMRI lists
%     The input fMRI data lists are assumed to exist in the following way:
%     For data in fsaverage space:
%     <out_dir>/data_list/fMRI_list/lh_sub?_sess?.txt
%     <out_dir>/data_list/fMRI_list/rh_sub?_sess?.txt
%     For data in fs_LR_32k space:
%     <out_dir>/data_list/fMRI_list/sub?_sess?.txt
%     Each line in above files corresponds to the full path of the fMRI
%     data of each run. For example, if the subject 1 in fsaverage5 has 2 
%     sessions and session 2 has 2 runs then there should be:
%     <out_dir>/data_list/fMRI_list/lh_sub1_sess1.txt
%     <out_dir>/data_list/fMRI_list/rh_sub1_sess1.txt
%     <out_dir>/data_list/fMRI_list/lh_sub1_sess2.txt
%     <out_dir>/data_list/fMRI_list/rh_sub1_sess2.txt
%     and lh_sub1_sess2.txt should have two lines:
%     <path_to_fMRI_data>/lh*<fMRI_filename_of_run_1>.nii.gz
%     <path_to_fMRI_data>/lh*<fMRI_filename_of_run_2>.nii.gz
%
%     2) censor lists
%     The input censor lists are assumed to exist in the following way:
%     <out_dir>/data_list/censor_list/sub?_sess?.txt
%     Each line in above file corresponds to the full path of the censor
%     list of each run. For example, if the subject 1 has 2 sessions and
%     session 2 has 2 runs then there should be:
%     <out_dir>/data_list/censor_list/sub1_sess1.txt
%     <out_dir>/data_list/censor_list/sub1_sess2.txt
%     and sub1_sess2.txt should have two lines:
%     <path_to_fMRI_data>/<censor_filename_of_run_1>
%     <path_to_fMRI_data>/<censor_filename_of_run_2>
%     Please note that the censor file shouled be a text file contains a
%     single column with binary numbers and its length is the number of 
%     timepoints. The outliers are indicated by 0s and will be flaged out
%     when compute the profiles. 
%
%     3) surface lists     
%     The input surface lists are assumed to be exist in the following way:
%     Some datasets such as HCP, might have the fs_LR_32k surface template 
%     in individual space:
%     ??????.L.midthickness.32k_fs_LR.surf.gii 
%     ??????.L.midthickness.32k_fs_LR.surf.gii 
%     If the user would like to pass in surface files, there should be
%     <out_dir>/data_list/surface_list/lh_sub?_sess?.txt;
%     Each line in this file corresponds to a single run of this subject.  
%     If there is no individual surface file, please leave the 
%     <out_dir>/data_list/surface_list as an exmpty folder.
%
%   - sub: (string)
%
%     The subject ID, e.g. '1','2'. 
%     The subject ID could be defined by 1,2,3... or other style 
%     0064,0021,0032... The user need to make sure the IDs are
%     consistent with the IDs they used in the fMRI lists and censor lists.
%
%   - num_sess: (string)
%
%     The total number of sessions, e.g. '4'. We use all sessions to generate
%     gradients.   
%
% Output:
%
%   The generated diffusion embedding matrices of gradients will saved into 
%   the following files:
%   <out_dir>/gradients/sub<sub_id>/lh_emb_100_distance_matrix.mat
%   <out_dir>/gradients/sub<sub_id>/rh_emb_100_distance_matrix.mat
%
% Example:
% CBIG_ArealMSHBM_generate_gradient('fsaverage6', '/path/project_dir/training_set','1','2');
% CBIG_ArealMSHBM_generate_gradient('fs_LR_32k', '/path/project_dir/training_set','1','4');
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
num_sess = str2double(num_sess);

out_gradient_dir = fullfile(out_dir,'gradients',['sub' sub]);
if(~exist(out_gradient_dir))
    mkdir(out_gradient_dir);
end
% check if there is any individual surface template
if(~exist(fullfile(out_dir,'data_list','surface_list'),'dir'))
    lh_ind_surf = 'NONE';
    rh_ind_surf = 'NONE';
    warning('Did not specify any individual surface template, will use the group average surface template.\n');
else
    dir_check = dir(fullfile(out_dir,'data_list','surface_list'));
    if(length(dir_check) == 2)%empty folder
        lh_ind_surf = 'NONE';
        rh_ind_surf = 'NONE';
        warning('Did not specify any individual surface template, will use the group average surface template.\n');
    else
        for sess = 1:num_sess
            lh_ind_surf_list = fullfile(out_dir,'data_list','surface_list',['lh_sub' sub '_sess' num2str(sess) '.txt']);
            rh_ind_surf_list = fullfile(out_dir,'data_list','surface_list',['rh_sub' sub '_sess' num2str(sess) '.txt']);
            curr_lh_surf = table2cell(readtable(lh_ind_surf_list, 'Delimiter',' ','ReadVariableNames',false));
            curr_rh_surf = table2cell(readtable(rh_ind_surf_list, 'Delimiter',' ','ReadVariableNames',false));
            if(sess == 1)
                lh_ind_surf = curr_lh_surf;
                rh_ind_surf = curr_rh_surf;
            else
                lh_ind_surf = [lh_ind_surf; curr_lh_surf];
                rh_ind_surf = [rh_ind_surf; curr_rh_surf];
            end
        end
    end
end
if(strcmp(mesh,'fs_LR_32k'))
    sub_FC = '10';
    sub_verts = '200';
    downsample = '3.2';
    medial_mask = 'NONE';

    for sess = 1:num_sess
        fMRI_list = fullfile(out_dir,'data_list','fMRI_list',['sub' sub '_sess' num2str(sess) '.txt']);
        censor_list = fullfile(out_dir,'data_list','censor_list',['sub' sub '_sess' num2str(sess) '.txt']);
        curr_fMRI = table2cell(readtable(fMRI_list, 'Delimiter',' ','ReadVariableNames',false));
        curr_censor = table2cell(readtable(censor_list, 'Delimiter',' ','ReadVariableNames',false));
        if(sess == 1)
            fMRI_all = curr_fMRI;
            censor_all = curr_censor;
        else
            fMRI_all = [fMRI_all; curr_fMRI];
            censor_all = [censor_all; curr_censor];
        end
    end
    if(~exist(fullfile(out_gradient_dir,['rh_gradient_distance_matrix.npy'])))
        if(~exist(fullfile(out_gradient_dir,['gradients_edge_density.dtseries.nii'])))
            CBIG_SPGrad_RSFC_gradients(fMRI_all, 'NONE', censor_all, lh_ind_surf, rh_ind_surf, ...
            mesh, medial_mask, sub_FC, sub_verts, out_gradient_dir);
        end
        CBIG_SPGrad_generate_gradient_matrix(mesh, medial_mask, downsample, out_gradient_dir);
    end
    cmd = ['sh ' CBIG_CODE_DIR '/utilities/matlab/speedup_gradients/CBIG_SPGrad_diffusion_embedding.sh '...
     out_gradient_dir ' 100' ];
    system(cmd); 
    CBIG_SPGrad_upsample_embed_matrix(mesh, medial_mask, 100, out_gradient_dir);            

elseif(~isempty(strfind(mesh,'fsaverage'))) 
    sub_FC = '100';
    sub_verts = '200';
    downsample = '3.2';
    medial_mask = 'NONE';
    for sess = 1:num_sess
        lh_fMRI_list = fullfile(out_dir,'data_list','fMRI_list',['lh_sub' sub '_sess' num2str(sess) '.txt']);
        rh_fMRI_list = fullfile(out_dir,'data_list','fMRI_list',['rh_sub' sub '_sess' num2str(sess) '.txt']);
        censor_list = fullfile(out_dir,'data_list','censor_list',['sub' sub '_sess' num2str(sess) '.txt']);
        curr_lh_fMRI = table2cell(readtable(lh_fMRI_list, 'Delimiter',' ','ReadVariableNames',false));
        curr_rh_fMRI = table2cell(readtable(rh_fMRI_list, 'Delimiter',' ','ReadVariableNames',false));
        curr_censor = table2cell(readtable(censor_list, 'Delimiter',' ','ReadVariableNames',false));
        if(sess == 1)
            lh_fMRI_all = curr_lh_fMRI;
            rh_fMRI_all = curr_rh_fMRI;
            censor_all = curr_censor;
        else
            lh_fMRI_all = [lh_fMRI_all; curr_lh_fMRI];
            rh_fMRI_all = [rh_fMRI_all; curr_rh_fMRI];
            censor_all = [censor_all; curr_censor];
        end
    end
    if(~exist(fullfile(out_gradient_dir,['rh_gradient_distance_matrix.npy'])))
        if(~exist(fullfile(out_gradient_dir,['gradients_edge_density.dtseries.nii'])))
            CBIG_SPGrad_RSFC_gradients(lh_fMRI_all, rh_fMRI_all, censor_all, lh_ind_surf, rh_ind_surf, ...
            mesh, medial_mask, sub_FC, sub_verts, out_gradient_dir);
        end
        CBIG_SPGrad_generate_gradient_matrix(mesh, medial_mask, downsample, out_gradient_dir);
    end
    cmd = ['sh ' CBIG_CODE_DIR '/utilities/matlab/speedup_gradients/CBIG_SPGrad_diffusion_embedding.sh '...
     out_gradient_dir ' 100' ];
    system(cmd); 
    CBIG_SPGrad_upsample_embed_matrix(mesh, medial_mask, 100, out_gradient_dir);   
else
    error('Unknown mesh type. Only fsaverage surface types, for example, fsaverage5/6 and fs_LR_32k will be allowed.')
end

