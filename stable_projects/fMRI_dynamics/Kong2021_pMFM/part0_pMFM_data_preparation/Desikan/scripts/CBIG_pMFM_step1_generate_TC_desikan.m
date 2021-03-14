function CBIG_pMFM_step1_generate_TC_desikan()

% This function is the wrapper to generate parcellated time serises and FCD
% for Desikan parcellation
%
% There is no input for this function as it can automatically get the
% output file from previous step.
% There is no output for this function as it will generate the output files
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

TC_result_dir = '../output/TC';
if ~exist(TC_result_dir,'dir')
    mkdir(TC_result_dir)
end

generate_training_TC()
generate_validation_TC()
generate_test_TC()

end


function generate_training_TC()

% This function is the wrapper to generate parcellated time serises and FCD
% for Desikan parcellation for training set
%
% There is no input for this function as it can automatically get the
% output file from previous step.
% There is no output for this function as it will generate the output files
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% loading the HCP directory
HCP_dir = getenv('CBIG_HCP_DIR');

%% loading training subject list
load('../input/SC_set.mat', 'sub_train')
train_list = sub_train;
run_list = {'rfMRI_REST1_LR','rfMRI_REST1_RL','rfMRI_REST2_LR','rfMRI_REST2_RL'};
TC_train_dir = '../output/TC/train';
if ~exist(TC_train_dir,'dir')
    mkdir(TC_train_dir)
end

%% loading DesikanKilliany36 labeling fslr32k
lh_mesh_fslr32k = CBIG_read_fslr_surface('lh','fs_LR_32k','inflated','aparc.annot');
rh_mesh_fslr32k = CBIG_read_fslr_surface('rh','fs_LR_32k','inflated','aparc.annot');

lh_fslr32k = lh_mesh_fslr32k.MARS_label;
rh_fslr32k = rh_mesh_fslr32k.MARS_label;

for k = 1:max(lh_fslr32k) %%% assuming continuous labeling
    lh_ROI_cell{k} = find(lh_fslr32k==k);
end
for k = 1:max(rh_fslr32k) %%% assuming continious labeling
    rh_ROI_cell{k} = find(rh_fslr32k==k);
end

%% TC for ROI testing
disp('Train part')
for i = 1:length(train_list)
    sub_index_train = num2str(train_list(i));
    for j = 1:4
        file_dir = [HCP_dir '/individuals/' sub_index_train '/MNINonLinear/Results/' run_list{j} '/'];
        file_name = [run_list{j} '_Atlas_MSMAll_hp2000_clean.dtseries.nii'];
        if exist([file_dir file_name],'file') == 0
            break;
        end
        cii_full = ft_read_cifti([file_dir file_name]);
        
        %% Extracting TC in lh
        lh_indexes = find(cii_full.brainstructure==1);
        lh_fullTC = cii_full.dtseries(lh_indexes,:);

        % Extract TC in Desikan's parcels
        lh_tc = zeros(length(lh_ROI_cell),size(lh_fullTC,2));
        for n = 1:length(lh_ROI_cell)
            lh_tc(n,:) = nanmean(lh_fullTC(lh_ROI_cell{n},:));
        end


        %% Extractiing TC in rh
        rh_indexes = find(cii_full.brainstructure==2);
        rh_fullTC = cii_full.dtseries(rh_indexes,:);

        % Extract TC in Desikan's parcels
        rh_tc = zeros(length(rh_ROI_cell),size(rh_fullTC,2));
        for n = 1:length(rh_ROI_cell)
            rh_tc(n,:) = nanmean(rh_fullTC(rh_ROI_cell{n},:));
        end
        
        TC = [lh_tc;rh_tc];
        save([TC_train_dir '/' sub_index_train '_' run_list{j} '.mat'],'TC')
    end
    disp(i)
end

end


function generate_validation_TC()

% This function is the wrapper to generate parcellated time serises and FCD
% for Desikan parcellation for validation set
%
% There is no input for this function as it can automatically get the
% output file from previous step.
% There is no output for this function as it will generate the output files
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% loading the HCP directory
HCP_dir = getenv('CBIG_HCP_DIR');

%% loading training subject list
load('../input/SC_set.mat', 'sub_vali')
vali_list = sub_vali;
run_list = {'rfMRI_REST1_LR','rfMRI_REST1_RL','rfMRI_REST2_LR','rfMRI_REST2_RL'};
TC_vali_dir = '../output/TC/validation';
if ~exist(TC_vali_dir,'dir')
    mkdir(TC_vali_dir)
end

%% loading DesikanKilliany36 labeling fslr32k
lh_mesh_fslr32k = CBIG_read_fslr_surface('lh','fs_LR_32k','inflated','aparc.annot');
rh_mesh_fslr32k = CBIG_read_fslr_surface('rh','fs_LR_32k','inflated','aparc.annot');

lh_fslr32k = lh_mesh_fslr32k.MARS_label;
rh_fslr32k = rh_mesh_fslr32k.MARS_label;

for k = 1:max(lh_fslr32k) %%% assuming continuous labeling
    lh_ROI_cell{k} = find(lh_fslr32k==k);
end
for k = 1:max(rh_fslr32k) %%% assuming continious labeling
    rh_ROI_cell{k} = find(rh_fslr32k==k);
end

%% TC for ROI testing
disp('Validation part')
for i = 1:length(vali_list)
    sub_index_vali = num2str(vali_list(i));
    for j = 1:4
        file_dir = [HCP_dir '/individuals/' sub_index_vali '/MNINonLinear/Results/' run_list{j} '/'];
        file_name = [run_list{j} '_Atlas_MSMAll_hp2000_clean.dtseries.nii'];
        if exist([file_dir file_name],'file') == 0
            break;
        end
        cii_full = ft_read_cifti([file_dir file_name]);
        
        %% Extracting TC in lh
        lh_indexes = find(cii_full.brainstructure==1);
        lh_fullTC = cii_full.dtseries(lh_indexes,:);

        % Extract TC in Desikan's parcels
        lh_tc = zeros(length(lh_ROI_cell),size(lh_fullTC,2));
        for n = 1:length(lh_ROI_cell)
            lh_tc(n,:) = nanmean(lh_fullTC(lh_ROI_cell{n},:));
        end


        %% Extractiing TC in rh
        rh_indexes = find(cii_full.brainstructure==2);
        rh_fullTC = cii_full.dtseries(rh_indexes,:);

        % Extract TC in Desikan's parcels
        rh_tc = zeros(length(rh_ROI_cell),size(rh_fullTC,2));
        for n = 1:length(rh_ROI_cell)
            rh_tc(n,:) = nanmean(rh_fullTC(rh_ROI_cell{n},:));
        end
        
        TC = [lh_tc;rh_tc];
        save([TC_vali_dir '/' sub_index_vali '_' run_list{j} '.mat'],'TC')
    end
    disp(i)
end

end


function generate_test_TC()

% This function is the wrapper to generate parcellated time serises and FCD
% for Desikan parcellation for test set
%
% There is no input for this function as it can automatically get the
% output file from previous step.
% There is no output for this function as it will generate the output files
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% loading the HCP directory
HCP_dir = getenv('CBIG_HCP_DIR');

%% loading training subject list
load('../input/SC_set.mat', 'sub_test')
test_list = sub_test;
run_list = {'rfMRI_REST1_LR','rfMRI_REST1_RL','rfMRI_REST2_LR','rfMRI_REST2_RL'};
TC_test_dir = '../output/TC/test';
if ~exist(TC_test_dir,'dir')
    mkdir(TC_test_dir)
end

%% loading DesikanKilliany36 labeling fslr32k
lh_mesh_fslr32k = CBIG_read_fslr_surface('lh','fs_LR_32k','inflated','aparc.annot');
rh_mesh_fslr32k = CBIG_read_fslr_surface('rh','fs_LR_32k','inflated','aparc.annot');

lh_fslr32k = lh_mesh_fslr32k.MARS_label;
rh_fslr32k = rh_mesh_fslr32k.MARS_label;

for k = 1:max(lh_fslr32k) %%% assuming continuous labeling
    lh_ROI_cell{k} = find(lh_fslr32k==k);
end
for k = 1:max(rh_fslr32k) %%% assuming continious labeling
    rh_ROI_cell{k} = find(rh_fslr32k==k);
end

%% TC for ROI testing
disp('test part')
for i = 1:length(test_list)
    sub_index_test = num2str(test_list(i));
    for j = 1:4
        file_dir = [HCP_dir '/individuals/' sub_index_test '/MNINonLinear/Results/' run_list{j} '/'];
        file_name = [run_list{j} '_Atlas_MSMAll_hp2000_clean.dtseries.nii'];
        if exist([file_dir file_name],'file') == 0
            break;
        end
        cii_full = ft_read_cifti([file_dir file_name]);
        
        %% Extracting TC in lh
        lh_indexes = find(cii_full.brainstructure==1);
        lh_fullTC = cii_full.dtseries(lh_indexes,:);

        % Extract TC in Desikan's parcels
        lh_tc = zeros(length(lh_ROI_cell),size(lh_fullTC,2));
        for n = 1:length(lh_ROI_cell)
            lh_tc(n,:) = nanmean(lh_fullTC(lh_ROI_cell{n},:));
        end


        %% Extractiing TC in rh
        rh_indexes = find(cii_full.brainstructure==2);
        rh_fullTC = cii_full.dtseries(rh_indexes,:);

        % Extract TC in Desikan's parcels
        rh_tc = zeros(length(rh_ROI_cell),size(rh_fullTC,2));
        for n = 1:length(rh_ROI_cell)
            rh_tc(n,:) = nanmean(rh_fullTC(rh_ROI_cell{n},:));
        end
        
        TC = [lh_tc;rh_tc];
        save([TC_test_dir '/' sub_index_test '_' run_list{j} '.mat'],'TC')
    end
    disp(i)
end

end


