function CBIG_preproc_compare_output(result_dir, ref_dir, output_dir)
% This function is used to compare the single subject preprocessing outputs
% from the entire pipeline. The volumes/surfaces at each stage will be 
% compared based on voxel time course correlation and absolute difference.
% This function is only used for unit-test purpose.
%
% Input:
%   -result_dir: 
%       the path for the preprocessing outputs generated from the
%       unit test
%   -ref_dir: 
%       the path for the refernce preprocessing outputs
%   -output_dir: 
%       the path to save out comparison results
%
% Ouput:
%   -CBIG_preproc_compare_output_result.txt:
%       a text file saved under 'output_dir' with the correlation and
%       maximum absolute difference between unit-test results and reference
%       results at each stage.
%
% [NOTE]: Please refer to prepro.config for the preprocessing steps used in
% this version of single subject preprocessing unit test.
% Written by Shaoshi Zhang under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'CBIG_fMRI_Preproc2016' ,'utilities'));
% create output file
fid = fopen(fullfile(output_dir, 'CBIG_preproc_compare_output_result.txt'),'wt');

% % loop through 2 runs (001, 002)
for i = [1, 2]
    run = ['00' num2str(i)];
    fprintf(fid, ['======Run ' run '======\n']);
    
    % raw image
    fprintf(fid, ['sub-NDARBF851NH6_bld' run '_rest\n']);
    result_path = fullfile(result_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest.nii.gz']);
    ref_path = fullfile(ref_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest.nii.gz']);
    compare_result(result_path, ref_path, fid);
    
    % skip8
    fprintf(fid, ['sub-NDARBF851NH6_bld' run '_rest_skip8\n']);
    result_path = fullfile(result_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest_skip8.nii.gz']);
    ref_path = fullfile(ref_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest_skip8.nii.gz']);    
    compare_result(result_path, ref_path, fid);

    
    % stc
    fprintf(fid, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc\n']);
    result_path = fullfile(result_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc.nii.gz']);
    ref_path = fullfile(ref_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc.nii.gz']);
    compare_result(result_path, ref_path, fid);
    
    % mc
    fprintf(fid, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc\n']);
    result_path = fullfile(result_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc.nii.gz']);
    ref_path = fullfile(ref_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc.nii.gz']); 
    compare_result(result_path, ref_path, fid);
    
    % sdc
    fprintf(fid, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc\n']);
    result_path = fullfile(result_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc.nii.gz']);
    ref_path = fullfile(ref_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc.nii.gz']); 
    compare_result(result_path, ref_path, fid);
    
    % regression with censor
    fprintf(fid, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc\n']);
    result_path = fullfile(result_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc.nii.gz']);
    ref_path = fullfile(ref_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc.nii.gz']); 
    compare_result(result_path, ref_path, fid);
    
    % interpolation
    fprintf(fid, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60\n']);
    result_path = fullfile(result_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60.nii.gz']);
    ref_path = fullfile(ref_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60.nii.gz']); 
    compare_result(result_path, ref_path, fid);
    
    % bandpass filtering
    fprintf(fid, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08\n']);
    result_path = fullfile(result_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08.nii.gz']);
    ref_path = fullfile(ref_dir, 'bold', run, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08.nii.gz']); 
    compare_result(result_path, ref_path, fid);
    
    % MNI2mm volume
    fprintf(fid, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_MNI2mm_sm6_finalmask\n']);
    result_path = fullfile(result_dir, 'vol', ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_MNI2mm_sm6_finalmask.nii.gz']);
    ref_path = fullfile(ref_dir, 'vol', ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_MNI2mm_sm6_finalmask.nii.gz']);
    compare_result(result_path, ref_path, fid);
    
    % fs5 surface (only correlation)
    fprintf(fid, ['sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_fs5\n']);    
    result_lh_path = fullfile(result_dir, 'surf', ['lh.sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_fs5']);
    result_rh_path = fullfile(result_dir, 'surf', ['rh.sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_fs5']);
    ref_lh_path = fullfile(ref_dir, 'surf', ['lh.sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_fs5']);
    ref_rh_path = fullfile(ref_dir, 'surf', ['rh.sub-NDARBF851NH6_bld' run '_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_fs5']);
    [~,~,corr_lh_ex,corr_rh_ex] = CBIG_preproc_compare_two_surfaces(result_lh_path, result_rh_path, ref_lh_path, ref_rh_path,'fsaverage5');
    fprintf(fid,'max correlation :%f\n',max([corr_lh_ex corr_rh_ex]));
    fprintf(fid,'min correlation :%f\n',min([corr_lh_ex corr_rh_ex]));
    fprintf(fid,'mean correlation :%f\n',mean([corr_lh_ex corr_rh_ex]));
    fprintf(fid,'median correlation :%f\n\n',median([corr_lh_ex corr_rh_ex]));
      
end

%FC matrix
fprintf(fid, ['====== FC ======\n']);

fprintf(fid, 'sub-NDARBF851NH6_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_all2all\n');
result_FC_path = fullfile(result_dir, 'FC_metrics', 'Pearson_r', 'sub-NDARBF851NH6_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_all2all.mat');
ref_FC_path = fullfile(ref_dir, 'FC_metrics', 'Pearson_r', 'sub-NDARBF851NH6_rest_skip8_stc_mc_sdc_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_all2all.mat');
result_FC_struct = load(result_FC_path);
ref_FC_struct = load(ref_FC_path);
result_FC = result_FC_struct.corr_mat;
ref_FC = ref_FC_struct.corr_mat;
result_FC_reshape = reshape(result_FC, size(result_FC, 1) * size(result_FC, 2), 1);
ref_FC_reshape = reshape(ref_FC, size(ref_FC, 1) * size(ref_FC, 2), 1);
FC_corr = corr(result_FC_reshape, ref_FC_reshape);
fprintf(fid,['FC correlation :' num2str(FC_corr) '\n']);
fprintf(fid,'FC max difference :%f\n\n',max(abs(result_FC(:) - ref_FC(:))));

fclose(fid);

    function vol_reshape = load_reshape(image_path)
    % This function loads an image and reshapes it to size(image,
    % 1)*size(image, 2) * size(image, 3) * frame

    % load an image and extract the volume
    vol_struct = MRIread(image_path);
    vol = vol_struct.vol;

    % reshape
    vol_reshape = reshape(vol, size(vol,1)*size(vol,2)*size(vol,3), size(vol,4));
    end


    function compare_result(result_path, ref_path, fid)
    % This function is used to calculate correlation and difference between
    % test image and reference image, and write the result to an output
    % file
    result_img = load_reshape(result_path);
    ref_img = load_reshape(ref_path);
    
    correlation = CBIG_preproc_corr_matrix(result_img, ref_img);
    fprintf(fid, ['max correlation: ' num2str(max(correlation(:))) '\n']);
    fprintf(fid, ['min correlation: ' num2str(min(correlation(:))) '\n']);
    fprintf(fid, ['mean correlation: ' num2str(mean(correlation(:))) '\n']);
    fprintf(fid, ['median correlation: ' num2str(median(correlation(:))) '\n']);
    fprintf(fid, ['max difference: ' num2str(max(abs(result_img(:) - ref_img(:)))) '\n\n']);              
    end


end