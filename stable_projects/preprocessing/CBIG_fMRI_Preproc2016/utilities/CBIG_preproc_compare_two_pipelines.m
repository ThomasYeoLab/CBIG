function CBIG_preproc_compare_two_pipelines(pipe_dir1, pipe_name1, pipe_stem1...
, pipe_dir2, pipe_name2, pipe_stem2, subject_id, run, output_dir, compare_type)

% CBIG_preproc_compare_two_pipelines(pipe_dir1, pipe_name1, pipe_stem1,...
% pipe_dir2, pipe_name2, pipe_stem2, subject_id, run, output_dir, compare_type)
%
% Given output dir, stem, subject id, run number, output_dir and compare_type of
% two preprocessing pipelines, this function will give you the correlation maps 
% and stats of these correlations between these two pipelines. For more details,
% read the input, output and example sections below.
%
% Input:
%   - pipe_dir1    :   
%     pipe_dir1 is the directory that contains all subjects' preprocessed folder
%     of pipeline1
%   - pipe_name1   :
%     name of pipeline1
%   - pipe_stem1   : 
%     pipeline1 output file stem. File stem is the string after run number but 
%     before the file extension. 
%     (e.g. '_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6_
%     fsaverage5' for surface output in procsurffast pipeline, 
%     '_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_FS1mm_FS2mm_sm6' for 
%     volume output in procsurffast pipeline)
%   - pipe_dir2    : 
%     pipe_dir2 is the directory that contains all subjects' preprocessed folder
%     of pipeline2
%   - pipe_name2   :
%     name of pipeline2
%   - pipe_stem2   : 
%     pipeline2 output file stem. File stem is the string after run number but 
%     before the file extension. 
%     (e.g. '_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6_
%     fsaverage5' for surface output in procsurffast pipeline, 
%     '_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_FS1mm_FS2mm_sm6' for 
%     volume output in procsurffast pipeline)
%   - subject_id   : subject name (e.g. 'Sub0068_Ses1')
%   - run          : bold run string. The bold run string must be three digit
%     number as shown in your [pipe_dir '/' subject_id '/bold'] (e.g. '002').
%     If you have multiple runs, please call current function multiple times.
%   - output_dir   : output directory
%   - compare_type : 'surf' OR 'vol'. If compare_type == 'surf', it's assumed 
%     that comparisons are done in fsaverage5 space; if compare_type == 'vol' 
%     && (pipe_name1 == 'procsurffast' || pipe_name2 == 'procsurffast'),
%     it's assumed that the comparisons are done in MNI 2mm space (128x128x128);
%     if compare_type == 'vol' && pipe_name1 ~= 'procsurffast' && pipe_name2 ~=
%     'procsurffast', it's assumed that comparisons are done in a small MNI 2mm
%     space (91x109x91).
%
% Output:
%   If compare_type == 'surf', the output file will be under the directory 
%   fullfile(output_dir, subject_id)
%   - [pipe_name1 '_' pipe_name2 '_corr_surf.png'] :
%     correlation map between surface data in pipe1 and pipe2
%   - [pipe_name1 '_' pipe_name2 '_corr_surf_hist.png'] :
%     histogram of correlation between surface data in pipe1 and pipe2
%   - [pipe_name1 '_' pipe_name2 '_corr_surf_stat.txt'] :
%     maximum, mean, minimum, median value of the correlation
%   If compare_type == 'vol', the output file will be under the directory
%   fullfile(output_dir, subject_id)
%   - [pipe_name1 '_' pipe_name2 '_corr_vol.nii.gz'] :
%     The correlation of all voxels in the whole brain between 2 pipelines are
%     saved in MNI 2mm volume space.
%   - [pipe_name1 '_' pipe_name2 '_corr_vol_gm.nii.gz'] :
%     The correlation of the voxels within grey matter mask between 2 pipelines
%     are saved in MNI2mm volume space.
%   - [pipe_name1 '_' pipe_name2 '_corr_vol_hist.png'] :
%     histogram of correlation of whole brain between 2 pipelines 
%   - [pipe_name1 '_' pipe_name2 '_corr_vol_gm_hist.png'] :
%     histogram of correlation of grey matter volume between 2 pipelines 
%   - [pipe_name1 '_' pipe_name2 '_corr_vol_stat.txt'] :
%     maximum, mean, minimum, median of correlation of whole brain between 
%     2 pipelines 
%   - [pipe_name1 '_' pipe_name2 '_corr_vol_gm_stat.txt'] :
%     maximum, mean, minimum, median of correlation of grey matter volume 
%     between 2 pipelines 
%
% Example:
% datadir1 = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
%    'CBIG_fMRI_Preproc2016', 'single_subject', 'data');
% datadir2 = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
%    'CBIG_fMRI_Preproc2016', '100subjects_clustering', 'preproc_out');
% CBIG_preproc_compare_two_pipelines(datadir1, 'pipename1', ...
% '_rest_reorient_skip_faln_mc_g1000000000_bpss_resid_fsaverage6_sm6_fsaverage5', ...
% datadir2, 'pipename2', '_rest_skip4_stc_mc_resid_lp0.08_fs6_sm6_fs5', 'Sub1116_Ses1', ...
% '002', fullfile(getenv('HOME'), 'storage', 'test'), 'surf')
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% create sub dir for each subject_id and each run
output_dir_sub_run = fullfile(output_dir, subject_id, run );
mkdir(output_dir_sub_run)

% get the full path of surface or volume in pipeline1
lh_surface_file1 = fullfile(pipe_dir1, subject_id, 'surf', ...
    ['lh.' subject_id '_bld' run pipe_stem1 '.nii.gz']);
rh_surface_file1 = fullfile(pipe_dir1, subject_id, 'surf', ...
    ['rh.' subject_id '_bld' run pipe_stem1 '.nii.gz']);
MNI2mm_vol_file1 = fullfile(pipe_dir1, subject_id, 'vol', ...
    [subject_id '_bld' run pipe_stem1 '.nii.gz']);

% get the full path of surface or volume in pipeline2
lh_surface_file2 = fullfile(pipe_dir2, subject_id, 'surf', ...
    ['lh.' subject_id '_bld' run pipe_stem2 '.nii.gz']);
rh_surface_file2 = fullfile(pipe_dir2, subject_id, 'surf', ...
    ['rh.' subject_id '_bld' run pipe_stem2 '.nii.gz']);
MNI2mm_vol_file2 = fullfile(pipe_dir2, subject_id, 'vol', ...
    [subject_id '_bld' run pipe_stem2 '.nii.gz']);

%%
% get MNI 2mm gm mask used to mask volume data
gm_mask_128_file = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'volume', ...
    'FSL_MNI152_masks', 'GM_Mask_MNI1mm_MNI2mm_128x128x128.nii.gz');
gm_mask_91_file = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'volume', ...
    'FSL_MNI152_masks/GM_Mask_MNI1mm_MNI2mm_91x109x91.nii.gz');

%% compare pipe_type1 and pipe_type2

outputname = [pipe_name1 '_' pipe_name2 '_corr'];

%%%%%%%%%%%%%%%%%%%%%
% compare the surface
%%%%%%%%%%%%%%%%%%%%%
if strcmp(compare_type,'surf') 
[corr_lh,corr_rh,corr_lh_ex,corr_rh_ex] = CBIG_preproc_compare_two_surfaces...
(lh_surface_file1,rh_surface_file1,lh_surface_file2,rh_surface_file2,'fsaverage5');

CBIG_DrawSurfaceMaps(corr_lh,corr_rh,'fsaverage5','inflated',0,1);
saveas(gcf,fullfile(output_dir_sub_run, [outputname '_surf.png']));
close(gcf);

figure
set(gcf,'Visible','off');
title('Surfaces correlation hist')
xlabel('Correlation between pipe1 and pipe2')
ylabel('Num of voxels')
hold on
hist([corr_lh_ex corr_rh_ex]);
hold off
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, fullfile(output_dir_sub_run, [outputname '_surf_hist.png']));

fid = fopen(fullfile(output_dir_sub_run, [outputname '_surf_stat.txt']),'w');
fprintf(fid,'max correlation :%f\n',max([corr_lh_ex corr_rh_ex]));
fprintf(fid,'mean correlation :%f\n',mean([corr_lh_ex corr_rh_ex]));
fprintf(fid,'min correlation :%f\n',min([corr_lh_ex corr_rh_ex]));
fprintf(fid,'median correlation :%f\n',median([corr_lh_ex corr_rh_ex]));
fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%
% compare MNI 2mm vol 
%%%%%%%%%%%%%%%%%%%%%
if strcmp(compare_type,'vol') 
[corr_vol,vol_size] = CBIG_preproc_compare_two_vols(MNI2mm_vol_file1,MNI2mm_vol_file2);

%load the gm mask and apply it to vol correlation
if strcmp(pipe_name1, 'procsurffast') || strcmp(pipe_name2, 'procsurffast')
    mri = MRIread(gm_mask_128_file);
else
    mri = MRIread(gm_mask_91_file);
end
mask=mri.vol;
mask1d=mask(:);
corr_vol_gm = corr_vol(logical(mask1d));

%write the vol correlation into nifti volume
corr_vol3d=reshape(corr_vol,size(mask));
mri.vol=corr_vol3d;
MRIwrite(mri,fullfile(output_dir_sub_run, [outputname '_vol.nii.gz']));

corr_gm_re=zeros(1,length(mask1d));
corr_gm_re(logical(mask1d))=corr_vol_gm;
corr_gm_vol=reshape(corr_gm_re,size(mask));
mri.vol=corr_gm_vol;
MRIwrite(mri,fullfile(output_dir_sub_run, [outputname '_vol_gm.nii.gz']));

%exclude the NaN value s 
corr_vol_ex = corr_vol(~isnan(corr_vol));
corr_vol_gm_ex = corr_vol_gm(~isnan(corr_vol_gm));

figure 
set(gcf,'Visible','off');
title('Volume correlation hist')
xlabel('Correlation between pipe1 and pipe2')
ylabel('Num of voxels')
hold on
hist(corr_vol_ex)
hold off
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,fullfile(output_dir_sub_run, [outputname '_vol_hist.png']));

figure 
set(gcf,'Visible','off');
title('Volume gm correlation hist')
xlabel('Correlation between pipe1 and pipe2')
ylabel('Num of voxels')
hold on
hist(corr_vol_gm_ex);
hold off
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,fullfile(output_dir_sub_run, [outputname '_vol_gm_hist.png']));

fid = fopen(fullfile(output_dir_sub_run, [outputname '_vol_stat.txt']),'w');
fprintf(fid,'max correlation :%f\n',max(corr_vol_ex));
fprintf(fid,'mean correlation :%f\n',mean(corr_vol_ex));
fprintf(fid,'min correlation :%f\n',min(corr_vol_ex));
fprintf(fid,'median correlation :%f\n',median(corr_vol_ex));
fclose(fid);

fid = fopen(fullfile(output_dir_sub_run, [outputname '_vol_gm_stat.txt']),'w');
fprintf(fid,'max correlation :%f\n',max(corr_vol_gm_ex));
fprintf(fid,'mean correlation :%f\n',mean(corr_vol_gm_ex));
fprintf(fid,'min correlation :%f\n',min(corr_vol_gm_ex));
fprintf(fid,'median correlation :%f\n',median(corr_vol_gm_ex));
fclose(fid);
end
