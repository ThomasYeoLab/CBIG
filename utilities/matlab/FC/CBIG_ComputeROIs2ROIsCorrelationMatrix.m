function CBIG_ComputeROIs2ROIsCorrelationMatrix(output_file, subj_text_list1, subj_text_list2, ...
    discard_frames_list, ROIs1, ROIs2, regression_mask1, regression_mask2, all_comb_bool, avg_sub_bool)

% CBIG_ComputeROIs2ROIsCorrelationMatrix(output_file, subj_text_list1, subj_text_list2, discard_frames_list, ...
% ROIs1, ROIs2, regression_mask1, regression_mask2, all_comb_bool, avg_sub_bool)
%   This function compute the ROIs to ROIs correlation matrix
%   The time courses for each subject will be grabbed and averaged within
%   the ROI. This averaged time course will be correlated with the time
%   course from another ROI.
%
%   -output_file: 
%   name of the output, it should end with '.mat'
%   e.g. output_file = 'lh2rh_corrmat.mat'
%
%   -subj_text_list1, subj_text_list2: 
%   text files in which each line
%   represents one subject's bold run
%   If a subject has multiple runs, it will still take one line, seperating
%   different runs by spaces 
%   
%   e.g.if subj_text_list1 contains 2 subjects and subject 2 has 2 runs,
%   then the file will look like: 
%   '<full_path>/subject1_run1_bold.nii.gz'
%   '<full_path>/subject2_run1_bold.nii.gz' <full_path>/subject2_run2_bold.nii.gz'
%   or
%   '<full_path>/subject1_run1_bold.dtseries.nii'
%   '<full_path>/subject2_run1_bold.dtseries.nii' <full_path>/subject2_run2_bold.dtseries.nii'

%   ROIs1 and ROIs2 can take a variable forms described below:
%
%   1. ROIs1 and ROIs2: text files in which each line corresponding to an ROI
%   Each ROI is a full path pointing to the location of that ROI (which
%   maybe a volume or a surface ROI)
%   
%   e.g., if ROIs1 contains 2 ROIs (.label files in this case) 
%   then the file will look like:
%   '<full_path>/lh.shen_region1.label'
%   '<full_path>/lh.shen_region2.label'
%   
%   2. ROIs1 and ROIs2 can also be a single .label/.nii.gz/.dlabel.nii/.annot file
%   It is useful when we want to compute the correlation of a single
%   surface ROI to a volume ROI
%   
%   e.g.
%   ROIs1 = '<full_path>/lh.shen_region1.label'
%   ROIs2 = '<full_path>/shen_3mm.nii.gz'
%   or
%   ROIs1 = '<full_path>/Thomas114_parcellation.dlabel.nii'
%   ROIs2 = '<full_path>/Thomas114_parcellation.dlabel.nii'
%   or
%   ROIs1 = '<full path>/lh.Schaefer2018_400Parcels_17Networks_order.annot'
%   ROIs2 = '<full path>/lh.Schaefer2018_400Parcels_17Networks_order.annot'
%   
%   -discard_frames_list: 
%   text file which has a structure like subj_text_list1 and subj_text_list2, 
%   i.e. each line will correspond to one subject but will point to the 
%   location of the frame index output of motion scrubing. The motion 
%   scrubing output should be a text file with a binary column, where 1 
%   means keep, 0 means throw.
%   
%   NOTE: discard_frames_list = 'NONE' if there are no motion scrubbing
%   
%   -regression_mask1, regression_mask2:
%   1. regression_mask1 and regression_mask2 can be a .mat file where each
%   of them contains a variable 'regress_cell' when we load this .mat file.
%   'regress_cell' is a 1xK cell. K is the number of regressors.
%   regress_cell{k} is the k-th regressor mask, which is a #voxels by 1 
%   binary vector.
%   
%   2. regression_mask1 and regression_mask2 can also be a 1xK cell, where
%   regress_cell{k} is the k-th regressor mask, which is a #voxels by 1 
%   binary vector.
%   
%   Signal within each mask is averaged to become a signal to be regressed 
%   from the fMRI volume/surface data.
%
%   If no regression is needed, then regress_list should be set to the string "NONE". 
%   The regression signals are jointly (rather than sequentially) regressed out.');
%
%   -all_comb_bool:
%   all_comb_bool = 1 if we want to compute all possible combination of
%   ROIs1 and ROIs2, that means if we have have N ROIs and M ROIs, our
%   output matrix will be M-by-N-by number of subjects
%
%   -all_sub_bool:
%   avg_sub_bool = 1 if we want our output matrix will be average across
%   all subjects
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(ischar(all_comb_bool))
    all_comb_bool = str2double(all_comb_bool);
end

if(ischar(avg_sub_bool))
    avg_sub_bool = str2double(avg_sub_bool);
end

%% read in subject_lists
subj_list_1 = read_sub_list(subj_text_list1);
subj_list_2 = read_sub_list(subj_text_list2);

if(length(subj_list_1) ~= length(subj_list_2))
   error('both lists should contain the same number of subjects'); 
end

%% read in discarded frames file name
if strcmp(discard_frames_list, 'NONE')
    frame = 0;
else
    frame = 1;
    discard_list = read_sub_list(discard_frames_list);
    if (length(discard_list) ~= length(subj_list_1))
        error('number of subjects in discard list should be the same as that of subjects');
    end
end 

%% read in ROIs
ROIs1_cell = read_ROI_list(ROIs1);
ROIs2_cell = read_ROI_list(ROIs2);

%% read in regression list
[regress_cell1, regress1] = read_regress_list(regression_mask1);
[regress_cell2, regress2] = read_regress_list(regression_mask2);

%% space allocation
if avg_sub_bool == 1
    if all_comb_bool == 1
        corr_mat = zeros(length(ROIs1_cell), length(ROIs2_cell));
    else
        if length(ROIs1_cell) ~= length(ROIs2_cell)
            error('ROIs1_cell should have the same length as ROIs2_cell when all_comb_bool = 0');
        end
        corr_mat = zeros(length(ROIs1_cell), 1);
    end % end if all_comb_bool = 1        
else
    if all_comb_bool == 1
        corr_mat = zeros(length(ROIs1_cell), length(ROIs2_cell), length(subj_list_1));
    else
        if length(ROIs1_cell) ~= length(ROIs2_cell)
            error('ROIs1_cell should have the same length as ROIs2_cell when all_comb_bool = 0');
        end
        corr_mat = zeros(length(ROIs1_cell), length(subj_list_1));
    end  
end % end if avg_sub_bool = 1

%% Compute correlation
for i = 1:length(subj_list_1) % loop through each subject
    disp(num2str(i));

    S1 = textscan(subj_list_1{i}, '%s');
    S1 = S1{1}; % a cell of size (#runs x 1) for subject i in the first list

    S2 = textscan(subj_list_2{i}, '%s');
    S2 = S2{1}; % a cell of size (#runs x 1) for subject i in the second list

    if frame
        discard = textscan(discard_list{i}, '%s');
        discard = discard{1}; % a cell of size (#runs x 1) for subject i in the discard list
    end

    for j = 1:length(S1)
        % retrieve the binary frame index
        if frame
            discard_file = discard{j};
            fid = fopen(discard_file, 'r');
            frame_index = fscanf(fid,'%d');
        end

        input = S1{j};
        if (isempty(strfind(input, '.dtseries.nii'))) % input is a nifti file: .nii.gz
            input_series = MRIread(input);
            % time_course1 will look like nframes x nvertices for e.g. 236 x 10242
            time_course1 = single(transpose(reshape(input_series.vol, ...
                size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), ...
                size(input_series.vol, 4))));
        else % input is a cifti file: .dtseries.nii
            input_series = ft_read_cifti(input);
            time_course1 = single(transpose(input_series.dtseries));
        end
        input = S2{j};
        if (isempty(strfind(input, '.dtseries.nii'))) % input is a nifti file: .nii.gz
            input_series = MRIread(input);
            % time_course1 will look like nframes x nvertices for e.g. 236 x 10242
            time_course2 = single(transpose(reshape(input_series.vol, ...
                size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), ...
                size(input_series.vol, 4))));
        else % input is a cifti file: .dtseries.nii
            input_series = ft_read_cifti(input);
            time_course2 = single(transpose(input_series.dtseries));
        end

        if frame
            time_course1(frame_index==0,:) = [];
            time_course2(frame_index==0,:) = [];
        end

        % create time_courses based on ROIs
        t_series1 = zeros(size(time_course1, 1), length(ROIs1_cell));
        for k = 1:length(ROIs1_cell)
            t_series1(:,k) = CBIG_nanmean(time_course1(:, ROIs1_cell{k}), 2);
        end

        t_series2 = zeros(size(time_course2, 1), length(ROIs2_cell));       
        for k = 1: length(ROIs2_cell)
            t_series2(:,k) = CBIG_nanmean(time_course2(:, ROIs2_cell{k}), 2);
        end

        % regression
        if(regress1)
            regress_signal = zeros(size(time_course1, 1), length(regress_cell1));
            for k = 1:length(regress_cell1)
               regress_signal(:, k) = CBIG_nanmean(time_course1(:, regress_cell1{k} == 1), 2); 
            end

            % faster than using glmfit in which we need to loop through
            % all voxels
            X = [ones(size(time_course1, 1), 1) regress_signal];
            pseudo_inverse = pinv(X);
            b = pseudo_inverse*t_series1;
            t_series1 = t_series1 - X*b;
        end

        if(regress2)
            regress_signal = zeros(size(time_course2, 1), length(regress_cell2));
            for k = 1:length(regress_cell2)
                regress_signal(:, k) = mean(time_course2(:, regress_cell2{k} == 1), 2);
            end

            % faster than using glmfit in which we need to loop through
            % all voxels
            X = [ones(size(time_course2, 1), 1) regress_signal];
            pseudo_inverse = pinv(X);
            b = pseudo_inverse*t_series2;
            t_series2 = t_series2 - X*b;
        end

        % normalize series (size of series now is nframes x nvertices)
        t_series1 = bsxfun(@minus, t_series1, mean(t_series1, 1));
        t_series1 = bsxfun(@times, t_series1, 1./sqrt(sum(t_series1.^2, 1)));

        t_series2 = bsxfun(@minus, t_series2, mean(t_series2, 1));
        t_series2 = bsxfun(@times, t_series2, 1./sqrt(sum(t_series2.^2, 1)));

        % compute correlation
        if all_comb_bool == 1
            subj_corr_mat = t_series1' * t_series2;
        else
            subj_corr_mat = transpose(sum(t_series1 .* t_series2, 1));
        end

        if j == 1
            subj_z_mat = CBIG_StableAtanh(subj_corr_mat); % Fisher r-to-z transform
        else
            subj_z_mat = subj_z_mat + CBIG_StableAtanh(subj_corr_mat);
        end            
    end % inner for loop for each run j of subject i

    subj_z_mat = subj_z_mat/length(S1); % average across number of runs

    if avg_sub_bool == 1
        corr_mat = corr_mat + subj_z_mat;
    else
        if all_comb_bool == 1
            corr_mat(:, :, i) = tanh(subj_z_mat);
        else
            corr_mat(:, i) = tanh(subj_z_mat);
        end
    end
end % outermost for loop
disp(['isnan: ' num2str(sum(isnan(corr_mat(:)))) ' out of ' num2str(numel(corr_mat))]);

if avg_sub_bool == 1
    corr_mat = corr_mat/length(subj_list_1);
    corr_mat = tanh(corr_mat);
end

%% write out results
if(~isempty(strfind(output_file, '.mat')))
    save(output_file, 'corr_mat', '-v7.3');
end
end

%% sub-function to read subject lists
function subj_list = read_sub_list(subject_text_list)
% this function will output a 1xN cell where N is the number of
% subjects in the text_list, each subject will be represented by one
% line in the text file
% NOTE: multiple runs of the same subject will still stay on the same
% line
% Each cell will contain the location of the subject, for e.g.
% '<full_path>/subject1_run1_bold.nii.gz <full_path>/subject1_run2_bold.nii.gz'
fid = fopen(subject_text_list, 'r');
i = 0;
while(1);
    tmp = fgetl(fid);
    if(tmp == -1)
        break
    else
        i = i + 1;
        subj_list{i} = tmp;
    end
end
fclose(fid);
end

%% sub-function to read ROI lists
function ROI_cell = read_ROI_list(ROI_list)
% ROI_list can be a .nii.gz/.mgz/.mgh/.dlabel.nii file contains a parcellation.

% this is for an arbitrary nii.gz file
if(~isempty(strfind(ROI_list, '.nii.gz')) || ~isempty(strfind(ROI_list, '.mgz')) || ~isempty(strfind(ROI_list, '.mgh')))
    ROI_vol = MRIread(ROI_list);

    regions = unique(ROI_vol.vol(ROI_vol.vol ~= 0));
    for i = 1:length(regions)
        ROI_cell{i} = find(ROI_vol.vol == regions(i));
    end
    
% this is for an arbitrary .dlabel.nii file
elseif (~isempty(strfind(ROI_list, '.dlabel.nii')))
    ROI_vol = ft_read_cifti(ROI_list, 'mapname','array');
    regions = unique(ROI_vol.dlabel(ROI_vol.dlabel ~= 0));
    for i = 1:length(regions)
        ROI_cell{i} = find(ROI_vol.dlabel == regions(i));
    end

% it can also be a single .label file
elseif (~isempty(strfind(ROI_list, '.label'))) % input ROI as a single .label file
    tmp = read_label([], ROI_list);
    ROI_cell{1} = tmp(:,1) + 1;
    
elseif (~isempty(strfind(ROI_list, '.annot'))) % input ROI as a single .annot file
    vertex_label = CBIG_read_annotation(ROI_list);
    regions = unique(vertex_label(vertex_label ~= 1));   % exclude medial wall
    for i = 1:length(regions)
        ROI_cell{i} = find(vertex_label == regions(i));
    end

else % input ROIs is a list of its locations: .nii.gz, dlabel.nii or .label
    fid = fopen(ROI_list, 'r');
    i = 0;
    while(1);
        tmp = fgetl(fid);
        if(tmp == -1)
            break
        else          
            if(~isempty(strfind(tmp, '.nii.gz')) || ~isempty(strfind(tmp, '.mgz')) || ~isempty(strfind(tmp, '.mgh')))
                ROI_vol = MRIread(tmp);
                regions = unique(ROI_vol.vol(ROI_vol.vol ~= 0));
                for n = 1:length(regions)
                    i = i + 1;
                    ROI_cell{i} = find(ROI_vol.vol == regions(n)); % each cell contains a list of vertex's indices
                end
            elseif(~isempty(strfind(tmp, '.dlabel.nii')))
                ROI_vol = ft_read_cifti(tmp, 'mapname','array');
                regions = unique(ROI_vol.dlabel(ROI_vol.dlabel ~= 0));
                for n = 1:length(regions)
                    i = i + 1;
                    ROI_cell{i} = find(ROI_vol.dlabel == regions(n));
                end
            elseif(~isempty(strfind(tmp, '.label')))
                i = i + 1;
                tmp = read_label([], tmp);
                ROI_cell{i} = tmp(:, 1) + 1;
            end
        end
    end % end while
    fclose(fid);
end

end

%% sub-function to read regression lists
function [regress_cell, isRegress] = read_regress_list(regression_mask)
if(strcmp(regression_mask, 'NONE'))
    isRegress = 0;
    regress_cell = [];
else
    isRegress = 1;
    %regression_cell is a existing variable
    if(exist('regression_cell', 'var') == 1)
        regress_cell = regression_mask;
    %regression_cell is a existing .mat file contains a variable regress_cell
    elseif (exist(regression_mask, 'file') == 2)
        load(regression_mask);
    else
        error('regression_cell should be a variable or a .mat file which contains a variable regress_cell');
    end

end

end
