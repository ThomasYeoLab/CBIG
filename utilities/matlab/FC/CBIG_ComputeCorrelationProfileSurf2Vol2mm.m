function CBIG_ComputeCorrelationProfileSurf2Vol2mm(seed_mesh, target, output_file, threshold, mask_file, input_file_txt, seed_file_txt, regress_bool, varargin)

% CBIG_ComputeCorrelationProfileSurf2Vol2mm(seed_mesh, target, output_file, threshold, mask_file, input_file_txt, seed_file_txt, regress_bool, varargin)
%
% seed_mesh:      resolution of seed regions
% target:         resolution of target regions
% output_file:    output file name
% threshold:      relative threshold of correlation matrix for binarization
% mask_file:      volumetric mask
% input_file_txt: a text file containing input volume data list
% seed_file_txt:  a text file containing seed surface data list
% regress_bool:   negative or 0 for no regression; 1 for Cerebellum and
%                 Striatum regression; 2 for Putamen regression
% varargin:       indices of subsets of ROIs (optional)
%
% Compute surface to volume correlation profiles with subcortical
% regression option. The regression masks are defined in FreeSurfer
% nonlinear space. Please see:
% (1) Buckner RL, Krienen FM, Castellanos A, Diaz JC, Yeo BTT. The
%     organization of the human cerebellum revealed by intrinsic functional
%     connectivity. J Neurophysiology, 106(5):2322-2345, 2011
% (2) Choi EY, Yeo BTT, Buckner RL. The organization of the human striatum
%     revealed by intrinsic functional connectivity. J Neurophysiology,
%     108(8):2242-2263, 2012  
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(nargin > 8)
   for i = 1:length(varargin)
        index(i) = str2num(varargin{i}); % only compute correlation for subsets of ROIs
   end
end

if(ischar(regress_bool))
   regress_bool = str2num(regress_bool); 
end

SUBJECTS_DIR = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'volume', 'FS_nonlinear_volumetric_space_4.5');
if(regress_bool > 0)
    
    % Load target voxels whose signals are to be cleaned
    if(regress_bool == 1)
        areas = {'Cerebellum', 'Striatum'}; 
        hemis = {'lh', 'rh'};
        cortex_dt = [3 4; ...
                     3.5 4.25];
        
        for k = 1:length(areas)
            
            if(strcmp(areas{k}, 'Cerebellum'))
                regress_target = MRIread([SUBJECTS_DIR '/LooseCerebellum.dist1.Mask.GCA.t0.5.nii.gz']);
            elseif(strcmp(areas{k}, 'Striatum'))
                regress_target = MRIread([SUBJECTS_DIR '/LooseStriatum.dist1.Mask.GCA.t0.5.nii.gz']);
            else
                error('Does not handle non cerebellum or striatum'); 
            end
            
            regress_target_cell{k} = find(regress_target.vol == 1); 
        end
    elseif(regress_bool == 2)
        areas = {'Putamen'}; 
        hemis = {'lh', 'rh'};
        cortex_dt = [4; ...
                     4.5];
                 
        regress_target = MRIread([SUBJECTS_DIR '/LooseStriatum.dist1.Mask.GCA.t0.5.nii.gz']);         
        regress_target_cell{1} = find(regress_target.vol == 1); 
    else
        error('No other regression for now!');
    end
    
    % load seeds voxels whose signals are to be regressed out
    for k = 1:length(areas)
        regress_seed_vol1 = MRIread([SUBJECTS_DIR '/scripts/SubcorticalRegression/CortexAdjacent.' hemis{1} '.' areas{k} '.GCA.DT' num2str(cortex_dt(1, k)) '.nii.gz']);
        regress_seed_vol2 = MRIread([SUBJECTS_DIR '/scripts/SubcorticalRegression/CortexAdjacent.' hemis{2} '.' areas{k} '.GCA.DT' num2str(cortex_dt(2, k)) '.nii.gz']);
        regress_seed = find(regress_seed_vol1.vol == 1 | regress_seed_vol2.vol == 1);
        [regress_seed(:, 1), regress_seed(:, 2), regress_seed(:, 3)] = ind2sub(size(regress_seed_vol1.vol), regress_seed);
        regress_seed_cell{k} = regress_seed;
    end
end

% read in seed files.
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

% read in input files.
fid = fopen(input_file_txt, 'r');
i = 0;
while(1);
   tmp = fscanf(fid, '%s\n', 1);
   if(isempty(tmp))
       break
   else
       i = i + 1;
       input_files{i} = tmp;
   end
end
fclose(fid);

mask = MRIread(mask_file);
mask_index = find(mask.vol == 1);

lh_seed_avg_mesh = CBIG_ReadNCAvgMesh('lh', seed_mesh, 'inflated', 'cortex');
rh_seed_avg_mesh = CBIG_ReadNCAvgMesh('rh', seed_mesh, 'inflated', 'cortex');

lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', target, 'inflated', 'cortex');
rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', target, 'inflated', 'cortex');

for i = 1:length(input_files)
    input_file = input_files{i};
    input_series = MRIread(input_file);
    t_series = transpose(reshape(input_series.vol, size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), size(input_series.vol, 4)));
    t_series = t_series(:, mask_index);
    
    % first perform regression if requested
    if(regress_bool > 0)
        disp(['Performing regression ' num2str(regress_bool)]);
        if(regress_bool == 1 || regress_bool == 2)
            for k = 1:length(areas)
                % extract signal from seed region (T x # cortex seeds)
                regress_seed = regress_seed_cell{k};
                regress_signal = zeros(size(input_series.vol, 4), size(regress_seed, 1));
                for jj = 1:size(regress_seed, 1)
                    regress_signal(:, jj) = input_series.vol(regress_seed(jj, 1), regress_seed(jj, 2), regress_seed(jj, 3), :);
                end

                % average cortex signal and regress entire subcortical structure
                regress_signal = mean(regress_signal, 2); % T x 1
                regress_target = regress_target_cell{k};

                [tmp, regress_target_index] = intersect(mask_index, regress_target);
                for kk = 1:length(regress_target_index)
                    regress_target_signal = t_series(:, regress_target_index(kk));
                    [b, dev, stats] = glmfit(regress_signal, regress_target_signal);
                    t_series(:, regress_target_index(kk)) = stats.resid;
                end
            end
        else
            error('No other regression for now!');
        end
    end
    
    seed_file = seed_files{i};
    hemi_index = strfind(seed_file, basename(seed_file));
    lh_seed_file = seed_file; lh_seed_file(hemi_index:hemi_index+1) = 'lh';
    lh_seed_series = MRIread(lh_seed_file);
    lh_seed_series = transpose(reshape(lh_seed_series.vol, size(lh_seed_series.vol, 1) * size(lh_seed_series.vol, 2) * size(lh_seed_series.vol, 3), size(lh_seed_series.vol, 4)));
    lh_seed_series = lh_seed_series(:, lh_avg_mesh.MARS_label(1:length(lh_seed_avg_mesh.MARS_label)) == 2);
    
    rh_seed_file = seed_file; rh_seed_file(hemi_index:hemi_index+1) = 'rh';
    rh_seed_series = MRIread(rh_seed_file);
    rh_seed_series = transpose(reshape(rh_seed_series.vol, size(rh_seed_series.vol, 1) * size(rh_seed_series.vol, 2) * size(rh_seed_series.vol, 3), size(rh_seed_series.vol, 4)));
    rh_seed_series = rh_seed_series(:, rh_avg_mesh.MARS_label(1:length(rh_seed_avg_mesh.MARS_label)) == 2);
    
    s_series = [lh_seed_series rh_seed_series];

    if(nargin > 8)
        s_series = s_series(:, index);
    end
    
    % normalize series (note that series are now of dimensions: T x N)
    disp('Computing Surf2Vol Profiles');
    s_series = bsxfun(@minus, s_series, mean(s_series, 1));
    s_series = bsxfun(@times, s_series, 1./sqrt(sum(s_series.^2, 1)));

    t_series = bsxfun(@minus, t_series, mean(t_series, 1));
    t_series = bsxfun(@times, t_series, 1./sqrt(sum(t_series.^2, 1)));
    
    corr_mat = s_series' * t_series;
    if(i == 1)
        output = corr_mat;
    else
        output = output + corr_mat;
    end
end
output = output / length(input_files);
corr_mat = output;
disp(['isnan: ' num2str(sum(isnan(corr_mat(:)))) ' out of ' num2str(numel(corr_mat))]);
corr_mat(isnan(corr_mat)) = 0;

tmp = sort(corr_mat(:), 'descend');
t = tmp(round(numel(corr_mat) * str2num(threshold)));
disp(['threshold: ' num2str(t)]);
if(str2num(threshold) < 1)
  corr_mat(corr_mat <  t) = 0;
  corr_mat(corr_mat >= t) = 1;
end

surf2vol_correlation_profile = transpose(corr_mat); % result is N voxels x # ROI 
save(output_file, 'surf2vol_correlation_profile', '-v7.3');

%exit
