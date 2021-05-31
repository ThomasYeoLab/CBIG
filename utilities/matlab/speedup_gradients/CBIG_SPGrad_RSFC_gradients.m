function CBIG_SPGrad_RSFC_gradients(fMRI_files1, fMRI_files2, censor_files, lh_ind_surf_file, ...
    rh_ind_surf_file, mesh, medial_mask, sub_FC, sub_verts, outputdir)

% CBIG_SPGrad_RSFC_gradients(fMRI_files1, fMRI_files2, censor_files, lh_ind_surf_file, ...
% rh_ind_surf_file, mesh, medial_mask, sub_FC, sub_verts, outputdir)
%
% This function generates resting-state function connectivity gradient maps for given fMRI data. 
% The gradinet algorithm was originally proposed by WashU group (Gordon et al. 2016). We implemented 
% a faster and less memory-intensive version by subsampling the functional connectivity matrices 
% (Kong et al. 2020). Furthermore, to reduce the memory usage, we split the RSFC similarity matrix
% into multiple small blocks (without doing any subsampling), where each block of the RSFC similarity
% matrix will be generated sequentially. Two subsampling parameters 'sub_FC' and 'sub_verts' control 
% how much data we will subsample.
%
% Input:
%     - fMRI_files1:
%       'fsaverage*' spaces:
%       the text file containing left hemisphere surface data list,
%       e.g. <path>/lh_fmri_data_list.txt;
%       Each line in this file corresponds to a single run of this subject. 
%       <fMRI_files1> can also be a path which points to a single run if the input
%       is fMRI data of a single run. 
%       e.g. <path>/lh.Sub0033_Ses1.nii.gz.
%       <fMRI_files1> can also be a 1x#run cell variable. Each element corresponds
%       to the path to a single run of this subject.
%
%       'fs_LR*' spaces:
%       the text file containing entire cortical surface data list,
%       e.g. <path>/fmri_data_list.txt;
%       Each line in this file corresponds to a single run of this subject.
%       <fMRI_files1> can also be a path which points to a single run if the input
%       is fMRI data of a single run.
%       e.g. <path>/Sub0033_Ses1.dtseries.nii.
%       <fMRI_files1> can also be a 1x#run cell variable. Each element corresponds
%       to the path to a single run of this subject.
%
%     - fMRI_files2:
%       'fsaverage*' spaces:
%       the text file containing right hemisphere surface data list,
%       e.g. <path>/rh_fmri_data_list.txt;
%       Each line in this file corresponds to a single run of this subject. 
%       <fMRI_files2> can also be a path which points to a single run if the input
%       is fMRI data of a single run. 
%       e.g. <path>/rh.Sub0033_Ses1.nii.gz.
%       <fMRI_files2> can also be a 1x#run cell variable. Each element corresponds
%       to the path to a single run of this subject.
%
%       'fs_LR*' spaces:
%       This input is not useful for input in 'fs_LR*' spaces (you can pass in 'NONE').
%
%     - lh_ind_surf_file, rh_ind_surf_file:
%       The text file containing path to the individual surface mesh file. 
%       Some datasets such as HCP, might have the fs_LR_32k surface template 
%       in individual space:
%       ??????.L.midthickness.32k_fs_LR.surf.gii 
%       ??????.L.midthickness.32k_fs_LR.surf.gii 
%       e.g. <path>/?h_surf_list.txt;
%       Each line in this file corresponds to a single run of this subject. 
%       <lh_ind_surf_file>/<rh_ind_surf_file> can also be a path which points to a
%       single run if the input is the surface mesh file of a single run. 
%       e.g. <path>/??????.?.midthickness.32k_fs_LR.surf.gii. 
%       <lh_ind_surf_file>/<rh_ind_surf_file> can also be a 1x#run cell variable.
%       Each element corresponds to the path to a single run of this subject.
%       If there is no individual surface file, the user can pass in 'NONE' for both
%       variables. 
%
%     - censor_files:
%       the text file containing outlier file name, 
%       e.g. <path>/outlier_censor_list.txt
%       Each line in this list corresponds to a single run of this subject.
%        <outlier_text> can also be a path which points to the outlier file of a single run 
%       if the input is fMRI data of a single run. 
%       e.g. <path>/outlier_Sub0033_Ses1_run1.txt
%       <censor_files> can also be a 1x#run cell variable. Each element corresponds
%       to the path to the outlier file of a single run of this subject.
%       Please note that the outlier file shouled be a text file contains a
%       single column with binary numbers and its length is the number of 
%       timepoints. The outliers are indicated by 0s and will be flaged out.
%
%     - mesh:
%       resolution of surface mesh, e.g. 'fsaverage6', 'fs_LR_32k'
%
%     - medial_mask: (#num_vertices x 1 binary vector or 'NONE')
%       the medial area mask. Set it as 'NONE' unless the medial wall area is defined differently
%       as fsaverage medial wall area from FreeSurfer or fs_LR_32k medial wall area from HCP. 
%       If the data has a different medial area mask, please pass in a #num_vertices x 1 binary
%       vector. The medial area should be denoted as 1, others should be denoted as 0.
%
%     - sub_FC, sub_verts:
%       Without any speed-up subsampling, the RSFC similarity matrix is generated by
%       correlation between two NxN RSFC matrix (N is number of vertices). The dimensionality
%       of the RSFC similarity matrix is therefore NxN. In the speed-up version, we subsample 
%       the RSFC matrices, so that RSFC similarity matrix is generated by correlation between 
%       a NxN2 RSFC matrix and a N2xN1 RSFC matrix. The dimensionality of the RSFC similarity
%       matrix is therefore NxN1. N2 = #vertices/<sub_FC>, N1 = #vertices/<sub_verts>.
%       Higher value indicates stronger subsampling.  
%       For fsaverage6, we suggest user set: sub_FC = '100'; sub_verts = '200'.
%       For fs_LR_32k, we suggest user set: sub_FC = '10'; sub_verts = '200'.
%
%     - output_dir:
%       output directory to save the results. The estimated gradient will be saved as
%       <output_dir>/gradients_edge_density.dtseries.nii
%
% Example:
%   lh_fMRI_files = '/data_list/fMRI_list/lh_sub1.txt';
%   rh_fMRI_files = '/data_list/fMRI_list/rh_sub1.txt';
%   censor_files = '/data_list/censor_list/sub1.txt';
%   lh_ind_surf = 'NONE';
%   rh_ind_surf = 'NONE';
%   sub_FC = '100';
%   sub_verts = '200';
%   medial_mask = 'NONE';
%   out_dir = '/Kong2020_ArealMSHBM/examples/generate_gradients/sub1';
%   CBIG_SPGrad_RSFC_gradients(lh_fMRI_fies, rh_fMRI_fies, censor_files, lh_ind_surf, rh_ind_surf, ...
%   'fsaverage6', medial_mask, sub_FC, sub_verts, out_dir);
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'),...
 '/external_packages/matlab/non_default_packages/cifti-matlab-WashU-gradient')));

%% SET UP ENVIROnMENT
disp('SET UP ENVIRONMENT ...')
% temporay output directory and dump file directory
outputtmpdir = [outputdir '/tmp'];
outputdumpdir = [outputdir '/dump'];

if (~exist(outputtmpdir))
    mkdir(outputtmpdir);
end
if (~exist(outputdumpdir))
    mkdir(outputdumpdir);
end

% midthickness surface and mesh structure 
if(strcmp(mesh,'fs_LR_32k'))
    atlasdir = fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'speedup_gradients',...
     'utilities', 'Conte69_atlas');
    lh_midsurf = fullfile(atlasdir, 'Conte69.L.midthickness.32k_fs_LR.surf.gii');
    rh_midsurf = fullfile(atlasdir, 'Conte69.R.midthickness.32k_fs_LR.surf.gii');
    % mesh structure
    lh_avg_mesh = CBIG_read_fslr_surface('lh', mesh, 'sphere', 'medialwall.annot');
    rh_avg_mesh = CBIG_read_fslr_surface('rh', mesh, 'sphere', 'medialwall.annot');

    if(strcmp(medial_mask,'NONE'))
        medial_mask = [lh_avg_mesh.MARS_label == 1;rh_avg_mesh.MARS_label == 1];
    end
    is_fsLR = 1;
elseif(~isempty(strfind(mesh,'fsaverage6')))
    atlasdir = fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'speedup_gradients',...
     'utilities', 'fs6_surface_template');
    lh_midsurf = fullfile(atlasdir, 'fsaverage6.L.midthickness.surf.gii');
    rh_midsurf = fullfile(atlasdir, 'fsaverage6.R.midthickness.surf.gii');

    % mesh structure
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'sphere', 'cortex');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'sphere', 'cortex');

    if(strcmp(medial_mask,'NONE'))
        medial_mask = [lh_avg_mesh.MARS_label==1 rh_avg_mesh.MARS_label==1]';
    end
    is_fsLR = 0;
else
    error(['We do not have midthickness surface mesh for ' mesh 'in current version.']);
end

%% SET UP PARAMETERS
% subsampling parameters for reducing memory and time requirement
sub_FC = str2num(sub_FC);
sub_verts = str2num(sub_verts);

% split the RSFC similarity matrix between FC_A and FC_B into multiple small blocks
% number of blocks for FC_A
a = 3;
% number of blocks for FC_B
b = 10;

% fisher z transform
dofisher = 0;

% sigma for geodesic smoothing applied to gradient maps
smooth = 2.55;

% K-vertices away neighbors for watershed algorithm local minima detection
K = 3;
% threshold parameters for watershed algorithm
minh = 100;
maxh = -100;

%% Read in fMRI and censor file list
% if fMRI file input is a text file contains fMRI file paths of multiple scans
if(isempty([strfind(fMRI_files1, '.dtseries.nii') strfind(fMRI_files1, '.nii.gz') strfind(fMRI_files1, '.mgh')]))
    fMRI_filenames1 = CBIG_SPGrad_read_sub_list(fMRI_files1);
    if(~is_fsLR)
        fMRI_filenames2 = CBIG_SPGrad_read_sub_list(fMRI_files2);
    end
    censor_filenames = CBIG_SPGrad_read_sub_list(censor_files);
    if(~strcmp(lh_ind_surf_file,'NONE'))
        lh_ind_surf = CBIG_SPGrad_read_sub_list(lh_ind_surf_file);
        if(length(lh_ind_surf) == 1)
            error(['lh_ind_surf_file only has one file,' ...
             'it should be a text file where each row is the path to the individual surface files.'])
        end
    else
        lh_ind_surf = repmat({lh_midsurf},1,length(fMRI_filenames1));
    end
    if(~strcmp(rh_ind_surf_file,'NONE'))
        rh_ind_surf = CBIG_SPGrad_read_sub_list(rh_ind_surf_file);
        if(length(rh_ind_surf) == 1)
            error(['rh_ind_surf_file only has one file,' ...
             'it should be a text file where each row is the path to the individual surface files.'])
        end
    else
        rh_ind_surf = repmat({rh_midsurf},1,length(fMRI_filenames1));
    end

% if fMRI file input is a cell structure
elseif(strcmp(class(fMRI_files1), 'cell'))
    fMRI_filenames1 = fMRI_files1;
    if(~is_fsLR)
        fMRI_filenames2 = fMRI_files2;
    end
    censor_filenames = censor_files;

    if(strcmp(lh_ind_surf_file,'NONE'))
        lh_ind_surf = repmat({lh_midsurf},1,length(fMRI_filenames1));
    else
        lh_ind_surf = lh_ind_surf_file;
    end
    if(strcmp(rh_ind_surf_file,'NONE'))
        rh_ind_surf = repmat({rh_midsurf},1,length(fMRI_filenames1));
    else
        rh_ind_surf = rh_ind_surf_file;
    end
% if fMRI file input is the fMRI file path of a single scan
else    
    fMRI_filenames1{1} = fMRI_files1;
    if(~is_fsLR)
        fMRI_filenames2{1} = fMRI_files2;
    end
    censor_filenames{1} = censor_files;

    if(strcmp(lh_ind_surf_file,'NONE'))
        lh_ind_surf{1} = lh_midsurf;
    else
        lh_ind_surf{1} = lh_ind_surf_file;
    end
    if(strcmp(rh_ind_surf_file,'NONE'))
        rh_ind_surf{1} = rh_midsurf;
    else
        rh_ind_surf{1} = rh_ind_surf_file;
    end
end

%% Neighbors for watershed algorithm local minima detection
rh_avg_mesh.vertexNbors(rh_avg_mesh.vertexNbors ~= 0) = ... 
rh_avg_mesh.vertexNbors(rh_avg_mesh.vertexNbors ~= 0) + size(lh_avg_mesh.vertexNbors, 2);

neighbors = double([lh_avg_mesh.vertexNbors rh_avg_mesh.vertexNbors]);
neighbors(neighbors == 0) = NaN;
neighbors = [1:size(neighbors, 2); neighbors];
neighbors = CBIG_SPGrad_neighbors_exclude_medial(neighbors, medial_mask);
neighbors = neighbors';
K_neighbors = CBIG_SPGrad_find_neighbors(neighbors, K);

%% GENERATE GRADIENTS
for scan_idx = 1:length(fMRI_filenames1)
    disp(['====================== Scan No.' num2str(scan_idx) '======================']);
    disp('Read in data ...')

    curr_fMRI1 = fMRI_filenames1{scan_idx};
    if(~is_fsLR)
        curr_fMRI2 = fMRI_filenames2{scan_idx};
    end
    curr_censor = censor_filenames{scan_idx};

    curr_tmask = load(curr_censor);
    [~, curr_data, ~] = CBIG_SPGrad_read_fmri(curr_fMRI1);
    curr_data(:, curr_tmask == 0) = [];
    
    if(sum(~isnan(curr_data(:))) == 0)
        error(['Current input data is not valid; all entries are NaN. Path: ' curr_fMRI1]);
    end
    curr_data(isnan(curr_data)) = 0;

    if(~is_fsLR)
        [~, rh_curr_data, ~] = CBIG_SPGrad_read_fmri(curr_fMRI2);
        rh_curr_data(:, curr_tmask == 0) = [];

        if(sum(~isnan(rh_curr_data(:))) == 0)
            error(['Current input data is not valid; all entries are NaN. Path: ' curr_fMRI2]);
        end
        rh_curr_data(isnan(rh_curr_data)) = 0;

        curr_data = [curr_data; rh_curr_data];
        clear rh_curr_data
    end

    % mask out medial wall and vertices with zero info
    curr_data(medial_mask, :) = [];
    num_vertices = size(curr_data, 1);

    % calculate number of vertices after downsampling
    disp('############################');
    disp('#Set up speed up parameters#');
    disp('############################');

    %setup random number generator
    rng(scan_idx, 'twister');

    [randinds_verts, randinds_FC] = CBIG_SPGrad_set_downsample_params(num_vertices, sub_verts, sub_FC);
    % num_sample_S is subsampled FC dimension N1 = #vertices/<sub_verts>
    num_sample_S = length(randinds_verts);
    % num_sample_FC is subsampled FC dimension N2 = #vertices/<sub_FC>
    num_sample_FC = length(randinds_FC);
    disp(['Downsample FC matrix FC_A with dimension: ' num2str(num_sample_S) ' x ' num2str(num_sample_FC)]);
    disp(['Downsample FC matrix FC_B with dimension: ' num2str(num_vertices) ' x ' num2str(num_sample_FC)]);

    disp('Split FC matrices FC_A and FC_B into small blocks ...')
    iter_a = a + (mod(num_sample_S, a) ~= 0);
    iter_b = b + (mod(num_vertices, b) ~= 0);
    block_size_a = fix(num_sample_S / a);
    block_size_b = fix(num_vertices / b);
    disp(['Split FC_A into ' num2str(iter_a) ' blocks, block size = ' num2str(block_size_a)]);
    disp(['Split FC_B into ' num2str(iter_b) ' blocks, block size = ' num2str(block_size_b)]);


    disp('####################################################');
    disp('#Compute FC similarity matrix between FC_A and FC_B#');
    disp('####################################################');
    %% calculate correlation FC matrix A (part of full correlation matrix)
    % subsample original time series with N2 random indices
    t_series = curr_data(randinds_FC, :);
    t_series = bsxfun(@minus, t_series, mean(t_series, 2));
    mag_t = sqrt(sum(t_series.^2, 2));

    for i = 1:iter_a
        disp(['============= FC_A: current block, iter_a ', num2str(i), '/', num2str(iter_a) '=============']);

        if (~exist(fullfile(outputtmpdir, num2str(i))))
            mkdir(fullfile(outputtmpdir, num2str(i)));
        end

        gradsname = ['gradients_LR_subcort_blockA_' num2str(i)];
        if (exist([outputdumpdir '/Done_grads_scan' num2str(scan_idx) '_block' num2str(i) '.dump']))
            disp('---------------------> Gradient file of this scan already generated for this block! Skip ...');
        else

            if (exist([outputtmpdir '/' num2str(i) '/FC_simi_LR_subcort.dtseries.nii']))
                disp('---------------------> FC simi file exist! Skip...');
            else
                % subsample original time series with N1 random indices, which were split into sequential blocks
                if (i == iter_a)
                    s_series = curr_data(randinds_verts(i * block_size_a - block_size_a + 1:end), :)';
                else
                    s_series = curr_data(randinds_verts(i * block_size_a - block_size_a + 1:i * block_size_a), :)';
                end

                s_series = bsxfun(@minus, s_series, mean(s_series, 1));
                mag_s = sqrt(sum(s_series.^2, 1));
                FC_A = (t_series * s_series) ./ (mag_t * mag_s);
                % apply the Fisher tranformation
                if (dofisher == 1)
                    FC_A = single(FisherTransform(FC_A));
                end

                clear s_series mag_s

                FC_simi_block = zeros(num_vertices, size(FC_A, 2), 'single');

                %% calculate correlation FC matrix B (part of full correlation matrix)
                for j = 1:iter_b
                    disp(['----> FC_B: current block, iter_b ', num2str(j), '/', num2str(iter_b)]);
                    tic;

                    if (j == iter_b)
                        s_series = curr_data(j * block_size_b - block_size_b + 1:end, :)';
                    else
                        s_series = curr_data(j * block_size_b - block_size_b + 1:j * block_size_b, :)';
                    end

                    s_series = bsxfun(@minus, s_series, mean(s_series, 1));
                    mag_s = sqrt(sum(s_series.^2, 1));
                    FC_B = (t_series * s_series) ./ (mag_t * mag_s);
                    % apply the Fisher tranformation
                    if (dofisher == 1)
                        FC_B = single(FisherTransform(FC_B));
                    end

                    clear s_series mag_s

                    %% compute FC similarity matrix between FC_A and FC_B
                    FC_A = bsxfun(@minus, FC_A, mean(FC_A, 1));
                    mag_a = sqrt(sum(FC_A.^2, 1));

                    FC_B = bsxfun(@minus, FC_B, mean(FC_B, 1));
                    mag_b = sqrt(sum(FC_B.^2, 1));

                    if (j == iter_b)
                        FC_simi_block(j * block_size_b - block_size_b + 1:end, :) = ...
                        single((FC_B' * FC_A) ./ (mag_b' * mag_a));
                    else
                        FC_simi_block(j * block_size_b - block_size_b + 1:j * block_size_b, :) = ...
                        single((FC_B' * FC_A) ./ (mag_b' * mag_a));
                    end

                    if (dofisher == 1)
                        FC_simi_block = FisherTransform(FC_simi_block);
                    end

                    clear FC_B mag_b

                    toc;
                end

                clear FC_A mag_a
                
                %clear FS_simi_block

                CBIG_SPGrad_construct_cifti_format(mesh, medial_mask, FC_simi_block, ...
                [outputtmpdir '/' num2str(i) '/FC_simi_LR_subcort.dtseries.nii']);

                clear FC_simi_block
            end

            disp('Next step ...')
            disp('###################################')
            disp('#Compute gradients using workbench#')
            disp('###################################')

            if (exist([outputtmpdir '/' num2str(i) '/' gradsname '.dtseries.nii']))
                disp('---------------------> Workbench gradient file exist! Skip...')
            else
                % calculate gradients
                system(['wb_command -cifti-gradient ' outputtmpdir '/' num2str(i) ...
                        '/FC_simi_LR_subcort.dtseries.nii COLUMN ' ...
                        outputtmpdir '/' num2str(i) '/' gradsname '.dtseries.nii -left-surface ' ...
                        lh_ind_surf{scan_idx} ' -right-surface ' rh_ind_surf{scan_idx}]);
            end

            % convert gradients and load
            grads = ft_read_cifti_mod([outputtmpdir '/' num2str(i) '/' gradsname '.dtseries.nii']);
            grads = single(grads.data);

            disp('Next step ...')
            disp('############################')
            disp('#Sum gradients across scans#')
            disp('############################')

            if (scan_idx == 1)
                grads_across_scan_block = grads;
                num_grads_across_scan_block = grads > 10^-10;
                save([outputtmpdir '/BLOCK' num2str(i) '_' gradsname '_across_scan.mat'], ...
                 'grads_across_scan_block', 'num_grads_across_scan_block', 'scan_idx','-v7.3');
                clear grads_across_scan_block num_grads_across_scan_block
            else
                tmp = load([outputtmpdir '/BLOCK' num2str(i) '_' gradsname '_across_scan.mat']);
                grads_across_scan_block = tmp.grads_across_scan_block + grads;
                num_grads_across_scan_block = tmp.num_grads_across_scan_block + (grads > 10^-10);
                save([outputtmpdir '/BLOCK' num2str(i) '_' gradsname '_across_scan.mat'], ...
                 'grads_across_scan_block', 'num_grads_across_scan_block', 'scan_idx','-v7.3');
                clear grads_across_scan_block num_grads_across_scan_block
            end

            % This dump file is used to check which scan is done. This is useful for resubmitting crashed jobs
            disp('Save dump file ...')
            CBIG_SPGrad_save_dump([outputdumpdir '/Done_grads_scan' num2str(scan_idx) '_block' num2str(i)]);
            rmdir([outputtmpdir '/' num2str(i)], 's');
        end

    end

end

disp('Next step ...')
rng(1, 'twister');
for i = 1:iter_a
    disp(['============= Gradients: current block, iter_a ', num2str(i), '/', num2str(iter_a) '=============']);

    gradsname = ['gradients_LR_subcort_blockA_' num2str(i)];
    avggradsname = ['gradients_LR_subcort_blockA_' num2str(i) '_avg_acrosssub'];

    disp('###############################')
    disp('#Average gradient across scans#')
    disp('###############################')

    if (exist([outputtmpdir '/BLOCK' num2str(i) '_' avggradsname '.dtseries.nii']))
        disp('---------------------> Average gradient file exist! Skip...');
    else
        load([outputtmpdir '/BLOCK' num2str(i) '_' gradsname '_across_scan.mat']);

        % average gradients across subjects
        grads_across_scan_block = grads_across_scan_block ./ num_grads_across_scan_block;
        grads_across_scan_block(num_grads_across_scan_block == 0) = 0; % in case the denominator is zero

        % save out average gradients
        CBIG_SPGrad_construct_cifti_format(mesh, medial_mask, grads_across_scan_block, ...
        [outputtmpdir '/BLOCK' num2str(i) '_' avggradsname '.dtseries.nii']);
        clear grads_across_scan_block num_grads_across_scan_block
        delete([outputtmpdir '/BLOCK' num2str(i) '_' gradsname '_across_scan.mat']);
    end

    % smooth gradients before edge detection
    disp('Next step ...')
    disp('############################')
    disp('#Smoothing average gradient#')
    disp('############################')

    if (exist([outputtmpdir '/BLOCK' num2str(i) '_' avggradsname '_smooth' num2str(smooth) '.dtseries.nii']))
        disp('---------------------> Smoothed gradient file exist! Skip..')
    else
        tic;
        system(['wb_command -cifti-smoothing ' ...
                outputtmpdir '/BLOCK' num2str(i) '_' avggradsname '.dtseries.nii ' num2str(smooth) ' ' ...
                num2str(smooth) ' COLUMN ' outputtmpdir '/BLOCK' num2str(i) '_' avggradsname '_smooth' ...
                num2str(smooth) '.dtseries.nii -left-surface ' lh_midsurf ...
                ' -right-surface ' rh_midsurf]);
        toc;
    end

    %% caculate edges
    disp('Next step ...')
    disp('#########################')
    disp('#Calculating water edges#')
    disp('#########################')

    if (exist([outputtmpdir '/BLOCK' num2str(i) '_' avggradsname '_smooth' num2str(smooth) '_edges.mat']))
        disp('---------------------> Watershed edge file exist! Skip...')
        load([outputtmpdir '/BLOCK' num2str(i) '_' avggradsname '_smooth' num2str(smooth) '_edges.mat']);
    else
        % load smoothed gradients
        disp('Compute minimametric ...')
        tic;
        tmp = ft_read_cifti_mod([outputtmpdir '/BLOCK' num2str(i) '_' avggradsname '_smooth' num2str(smooth)...
         '.dtseries.nii']);
        sm_avg_grads_across_scan_block = tmp.data;
        clear tmp

        if (min(sm_avg_grads_across_scan_block(:)) <= minh)
            minh = min(sm_avg_grads_across_scan_block(:));
        end

        if (max(sm_avg_grads_across_scan_block(:)) >= maxh)
            maxh = max(sm_avg_grads_across_scan_block(:));
        end

        % get local minima of each smoothed gradient map
        minimametrics = CBIG_SPGrad_find_minima(sm_avg_grads_across_scan_block, K_neighbors);
        toc;
        disp('Compute edges ...')
        tic;
        edges = watershed_algorithm_all_par_cifti(sm_avg_grads_across_scan_block, minimametrics, 50, 1,...
         neighbors, minh, maxh);
        toc;
        save([outputtmpdir '/BLOCK' num2str(i) '_' avggradsname '_smooth' num2str(smooth) '_edges.mat'], ...
         'edges','minimametrics')

        clear sm_avg_grads_across_scan_block minimametrics
    end

    if (i == 1)
        sum_edges = sum(edges == 0, 2);
    else
        sum_edges = sum_edges + sum(edges == 0, 2);
    end

end

% average across gradient maps and save
edge_density = sum_edges / num_sample_S;
CBIG_SPGrad_construct_cifti_format(mesh, medial_mask, edge_density, [outputdir '/gradients_edge_density.dtseries.nii']);
disp('########## Done!')

rmpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), ...
'/external_packages/matlab/non_default_packages/cifti-matlab-WashU-gradient')));