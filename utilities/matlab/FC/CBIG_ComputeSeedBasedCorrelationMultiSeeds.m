function output = CBIG_ComputeSeedBasedCorrelationMultiSeeds(seed_file_txt, target_file_txt, seed_mask_txt, target_mask_file, output_file)

% output = CBIG_ComputeSeedBasedCorrelationMultiSeeds(seed_file_txt, target_file_txt, seed_mask_txt, target_mask_file, output_file)
%
% seed_file_txt is list of files that are time series to be used as seed regions
% target_file_txt is list of files that are time series to be used as target regions
% seed_mask_txt is a list of files, each of which is a binary mask defining the seed
% target_mask_file is a single file that is a binary mask defining the region of target time series to perform corrrelation
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



% Read in target and seed files
if(strfind(seed_file_txt, '.txt'))
    fid = fopen(seed_file_txt, 'r');
    i = 0;
    while(1);
        tmp = fscanf(fid, '%s\n', 1);
        if(isempty(tmp))
            break
        else
            i = i + 1;
            seed_files{i} = tmp;
        end
    end
    fclose(fid);
else
    seed_files{1} = seed_file_txt;
end

if(strfind(target_file_txt, '.txt'))
    fid = fopen(target_file_txt, 'r');
    i = 0;
    while(1);
        tmp = fscanf(fid, '%s\n', 1);
        if(isempty(tmp))
            break
        else
            i = i + 1;
            target_files{i} = tmp;
        end
    end
    fclose(fid);
else
    target_files{1} = target_file_txt;
end

% read in seed index files
if(strfind(seed_mask_txt, '.txt'))
    fid = fopen(seed_mask_txt, 'r');
    i = 0;
    while(1);
        tmp = fscanf(fid, '%s\n', 1);
        if(isempty(tmp))
            break
        else
            i = i + 1;
            seed_mask_files{i} = tmp;
        end
    end
    fclose(fid);
else
    seed_mask_files{1} = seed_mask_txt;
end

for i = 1:length(seed_mask_files)
    if(strfind(seed_mask_files{i}, '.label'))
        label = read_label([], seed_mask_files{i});
        seed_indices{i} = label(:, 1) + 1;
    else
        seed_mask = MRIread(seed_mask_files{i});
        seed_indices{i} = find(seed_mask.vol == 1);
    end
end

% read target mask
if(strfind(target_mask_file, '.label'))
    label = read_label([], target_mask_file);
    target_index = label(:, 1) + 1;
else
    target_mask = MRIread(target_mask_file);
    target_index = find(target_mask.vol == 1);
end

disp('Seed Based Correlation');
tic
for i = 1:length(seed_files)
   
    disp(num2str(i));
    % read in time series
    seed_file = seed_files{i};
    seed_series = MRIread(seed_file);
    orig_s_series = transpose(reshape(seed_series.vol, size(seed_series.vol, 1) * size(seed_series.vol, 2) * size(seed_series.vol, 3), size(seed_series.vol, 4)));
    
    target_file = target_files{i};
    target_series = MRIread(target_file);
    t_series = transpose(reshape(target_series.vol, size(target_series.vol, 1) * size(target_series.vol, 2) * size(target_series.vol, 3), size(target_series.vol, 4)));
    t_series = t_series(:, target_index);
    
    % normalize series (note that series are now of dimensions: T x N)
    t_series = t_series - repmat(mean(t_series, 1), size(t_series, 1), 1);
    t_series = t_series./repmat(sqrt(sum(t_series.^2, 1)), size(t_series, 1), 1);
    
    % create empty frame
    if(i == 1)
        frame = zeros(size(target_series.vol, 1), size(target_series.vol, 2), size(target_series.vol, 3));
    end
    
    % Perform correlation for each seed
    for j = 1:length(seed_indices)
    
        % extract seed signal
        s_series = mean(orig_s_series(:, seed_indices{j}), 2);
        s_series = s_series - repmat(mean(s_series, 1), size(s_series, 1), 1);
        s_series = s_series./repmat(sqrt(sum(s_series.^2, 1)), size(s_series, 1), 1);
       
        % Compute correlation
        corr_val = sum(repmat(s_series, 1, size(t_series, 2)) .* t_series, 1);
        frame(target_index) = transpose(corr_val);
        if(i == 1 && j == 1)
            output = target_series;
            output.vol = zeros(size(target_series.vol, 1), size(target_series.vol, 2), size(target_series.vol, 3), length(seed_indices));
            output.vol(:, :, :, j) = frame; 
        else
            output.vol(:, :, :, j) = output.vol(:, :, :, j) + frame;
        end
    end
end
output.vol = output.vol/length(seed_files);
toc

% Write out results
if(nargin >= 5)
    MRIwrite(output, output_file);
    exit;
end








