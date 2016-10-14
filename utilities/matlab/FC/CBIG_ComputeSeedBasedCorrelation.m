function output = CBIG_ComputeSeedBasedCorrelation(seed_file_txt, target_file_txt, seed_mask_file, target_mask_file, type, output_file)

% output = CBIG_ComputeSeedBasedCorrelation(seed_file_txt, target_file_txt, seed_mask_file, target_mask_file, type, output_file)
%
% seed_file_txt is list of files that are time series to be used as seed regions
% target_file_txt is list of files that are time series to be used as target regions
% seed_mask_file is a single file that is a binary mask defining the seed
% target_mask_file is a single file that is a binary mask defining the region of target time series to perform corrrelation
%
% type can either be "average" or "full": "average" means average results
% across files. "full" means keep the full results.
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

% read in index files
if(strfind(seed_mask_file, '.label'))
    label = read_label([], seed_mask_file);
    seed_index = label(:, 1) + 1;
else
    seed_mask = MRIread(seed_mask_file);
    seed_index = find(seed_mask.vol == 1);
end

if(strfind(target_mask_file, '.label'))
    label = read_label([], target_mask_file);
    target_index = label(:, 1) + 1;
else
    target_mask = MRIread(target_mask_file);
    target_index = find(target_mask.vol == 1);
end

disp('Seed Based Correlation');
for i = 1:length(seed_files)
   
    % read in time series
    seed_file = seed_files{i};
    seed_series = MRIread(seed_file);
    s_series = transpose(reshape(seed_series.vol, size(seed_series.vol, 1) * size(seed_series.vol, 2) * size(seed_series.vol, 3), size(seed_series.vol, 4)));
    s_series = mean(s_series(:, seed_index), 2);
    
    target_file = target_files{i};
    target_series = MRIread(target_file);
    t_series = transpose(reshape(target_series.vol, size(target_series.vol, 1) * size(target_series.vol, 2) * size(target_series.vol, 3), size(target_series.vol, 4)));
    t_series = t_series(:, target_index);
    
    % normalize series (note that series are now of dimensions: T x N)
    s_series = bsxfun(@minus, s_series, mean(s_series, 1));
    s_series = bsxfun(@times, s_series, 1./sqrt(sum(s_series.^2, 1)));
    
    t_series = bsxfun(@minus, t_series, mean(t_series, 1));
    t_series = bsxfun(@times, t_series, 1./sqrt(sum(t_series.^2, 1)));
    
    corr_val = sum(bsxfun(@times, t_series, s_series), 1);
    
    if(i == 1)
        output = target_series;
        if(strcmp(type, 'average'))
            output.vol = zeros(size(target_series.vol, 1), size(target_series.vol, 2), size(target_series.vol, 3));
            output.vol(target_index) = transpose(corr_val);
        elseif(strcmp(type, 'full'))
            frame = zeros(size(target_series.vol, 1), size(target_series.vol, 2), size(target_series.vol, 3));
            output.vol = zeros(size(target_series.vol, 1), size(target_series.vol, 2), size(target_series.vol, 3), length(seed_files));
            frame(target_index) = transpose(corr_val);
            output.vol(:,:,:, i) = frame;
        else
            error(['Type ' type ' not recognized']);
        end
    else
        if(strcmp(type, 'average'))
            output.vol(target_index) = output.vol(target_index) + transpose(corr_val);
        elseif(strcmp(type, 'full'))
            frame(target_index) = transpose(corr_val);
            output.vol(:,:,:, i) = frame;
        else
            error(['Type ' type ' not recognized']); 
        end
    end
end

if(strcmp(type, 'average'))
    output.vol = output.vol/length(seed_files);
end
    
% Write out results
if(nargin >= 6)
    MRIwrite(output, output_file);
    exit;
end








