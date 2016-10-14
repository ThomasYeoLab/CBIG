function CBIG_ComputeFullVol2SurfCorrelationIndex(output_dir, lh_varargin_text, rh_varargin_text, vol_varargin_text, seed_mask_file, regress_bool, start_index, stop_index, pval)

% CBIG_ComputeFullVol2SurfCorrelationIndex(output_dir, lh_varargin_text, rh_varargin_text, vol_varargin_text, seed_mask_file, regress_bool, start_index, stop_index, pval)
%
% output_file:    output file name
% seed_mask_file:      surface mask
% lh_varargin_text, rh_varargin_text: a text file containing input surface 
%                                     data list of left and right
%                                     hemisphere.
% vol_varargin_text:  a text file containing volumetric data list
% regress_bool:   negative or 0 for no regression; 1 for Cerebellum and
%                 Striatum regression; 2 for Putamen regression
% start_index, stop_index: only consider indices from start_index to 
%                          stop_index that have not been already produced
% pval: If pval > 0 then this function will compute statistics
%
% Compute full surface to volume correlation profiles with subcortical
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


% read in both left and right text files.
fid = fopen(lh_varargin_text, 'r');
i = 0;
while(1);
   tmp = fgetl(fid);
   if(tmp == -1)
       break
   else
       i = i + 1;
       lh_varargin{i} = tmp;
   end
end
fclose(fid);

% read in both left and right text files.
fid = fopen(rh_varargin_text, 'r');
i = 0;
while(1);
   tmp = fgetl(fid);
   if(tmp == -1)
       break
   else
       i = i + 1;
       rh_varargin{i} = tmp;
   end
end
fclose(fid);

% read in volume
fid = fopen(vol_varargin_text, 'r');
i = 0;
while(1);
   tmp = fgetl(fid);
   if(tmp == -1)
       break
   else
       i = i + 1;
       vol_varargin{i} = tmp;
   end
end
fclose(fid);

% read seed_mask
seed_mask = MRIread(seed_mask_file);
seed_mask_index = find(seed_mask.vol(:) == 1);

% only consider indices from start_index to stop_index that have not been already produced
if(ischar(start_index))
   start_index = str2num(start_index);
end

if(ischar(stop_index))
    stop_index = str2num(stop_index);
end

if(nargin < 9)
    pval = 0;
else
    if(ischar(pval))
        pval = str2num(pval);
    end
end

seed_mask.vol(:) = 0;
for i = min(start_index, length(seed_mask_index)):min(stop_index, length(seed_mask_index))

    current_index = seed_mask_index(i);
    [y,x,z] = ind2sub(size(seed_mask.vol), current_index);
    y = y - 1;
    x = x - 1;
    z = z - 1;

    if(~exist(fullfile(output_dir, ['rh.' num2str(x) '_' num2str(y) '_' num2str(z) '.corr.nii.gz']), 'file'))
        seed_mask.vol(current_index) = 1;
    end
end
seed_mask_index = find(seed_mask.vol == 1);

    log_file = [vol_varargin_text '.log'];
    delete(log_file);
    system(['echo num_seeds: ' num2str(length(seed_mask_index)) ' >> ' log_file]);
    disp(['num_seeds: ' num2str(length(seed_mask_index))]); 

if(~isempty(seed_mask_index))

    % handle regression
    if(ischar(regress_bool))
        regress_bool = str2num(regress_bool);
    end

    SUBJECTS_DIR = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'volume', 'FS_nonlinear_volumetric_space_4.5');
    if(regress_bool > 0)

        areas = {'Cerebellum', 'Striatum'};
        hemis = {'lh', 'rh'};
        cortex_dt = [3 4; ...
            3.5 4.25];

        % Load target voxels whose signals are to be cleaned
        if(regress_bool == 1)
            for k = 1:length(areas)

                if(strcmp(areas{k}, 'Cerebellum'))
                    regress_target = MRIread([SUBJECTS_DIR '/scripts/downsample_data/LooseCerebellum.dist1.Mask.GCA.t0.5.nii.gz']);
                elseif(strcmp(areas{k}, 'Striatum'))
                    regress_target = MRIread([SUBJECTS_DIR '/scripts/downsample_data/LooseStriatum.dist1.Mask.GCA.t0.5.nii.gz']);
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

            regress_target = MRIread([SUBJECTS_DIR '/scripts/downsample_data/LooseStriatum.dist1.Mask.GCA.t0.5.nii.gz']);
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

    % Compute correlation
    for i = 1:length(lh_varargin)
        
        system(['echo ' num2str(i) ' >> ' log_file]);
        disp(num2str(i));

        C1 = textscan(lh_varargin{i}, '%s');
        C1 = C1{1};

        C2 = textscan(rh_varargin{i}, '%s');
        C2 = C2{1};

        C3 = textscan(vol_varargin{i}, '%s');
        C3 = C3{1};

        for j = 1:length(C1)

            % read surface time series
            input = C1{j};
            input_series = MRIread(input);
            t_series1 = single(transpose(reshape(input_series.vol, size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), size(input_series.vol, 4))));
            
            if(i == 1 && j == 1)
               output_mask = input_series; 
               output_mask.vol = zeros(size(output_mask.vol, 1), size(output_mask.vol, 2), size(output_mask.vol, 3));
               output_mask.nframes = 1;
               lh_num_vertices = numel(output_mask.vol); 
            end
            
            input = C2{j};
            input_series = MRIread(input);
            t_series2 = single(transpose(reshape(input_series.vol, size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), size(input_series.vol, 4))));

            t_series = [t_series1 t_series2];
            
            % read volume time series
            input = C3{j};
            input_series = MRIread(input);
            s_series = single(transpose(reshape(input_series.vol, size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), size(input_series.vol, 4))));
            s_series = s_series(:, seed_mask_index);

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

                        [tmp, regress_target_index] = intersect(seed_mask_index, regress_target);
                        for kk = 1:length(regress_target_index)
                            regress_target_signal = s_series(:, regress_target_index(kk));
                            [b, dev, stats] = glmfit(regress_signal, regress_target_signal);
                            s_series(:, regress_target_index(kk)) = stats.resid;
                        end
                    end
                else
                    error('No other regression for now!');
                end
            end

            % normalize series (note that series are now of dimensions: T x N)
            t_series = bsxfun(@minus, t_series, mean(t_series, 1));
            t_series = bsxfun(@times, t_series, 1./sqrt(sum(t_series.^2, 1)));

            s_series = bsxfun(@minus, s_series, mean(s_series, 1));
            s_series = bsxfun(@times, s_series, 1./sqrt(sum(s_series.^2, 1)));

            corr_mat = t_series' * s_series;
            if(j == 1)
                sbj_z_mat = CBIG_StableAtanh(corr_mat); % fisher-z transform
            else
                sbj_z_mat = sbj_z_mat + CBIG_StableAtanh(corr_mat);
            end
        end
        sbj_z_mat = sbj_z_mat/length(C1);
        clear corr_mat;

        disp(['isnan: ' num2str(sum(isnan(sbj_z_mat(:)))) ' out of ' num2str(numel(sbj_z_mat))]);
        sbj_z_mat(isnan(sbj_z_mat)) = 0;
        
        if(i == 1)
            output = sbj_z_mat;
            if(pval > 0)
                zmat = sbj_z_mat;
                zsqmat = sbj_z_mat.^2;
            end
        else
            output = output + sbj_z_mat;
            if(pval > 0)
                zmat = zmat + sbj_z_mat;
                zsqmat = zsqmat + sbj_z_mat.^2;
            end
        end
        
        clear sbj_z_mat;
    end
    output = output / length(lh_varargin);
    corr_mat = tanh(output);
    disp(['isnan: ' num2str(sum(isnan(corr_mat(:)))) ' out of ' num2str(numel(corr_mat))]);
    corr_mat(isnan(corr_mat)) = 0;
    clear output;

    if(pval > 0)
        system(['echo Computing statistics >> ' log_file]);
        disp('Computing statistics');

        zmat = zmat/length(lh_varargin); % compute mean
        zsqmat = zsqmat/length(lh_varargin);
        stdz = zsqmat - zmat.^2; stdz(stdz < 0) = 0;
        stdz = sqrt(stdz)*length(lh_varargin)/(length(lh_varargin) - 1); %compute std
        clear zsqmat;
        tmat = sqrt(length(lh_varargin))*zmat./(stdz + eps); % compute t stats
        clear stdz; clear zmat;
        pmat = -log10(tcdf(-abs(tmat), length(lh_varargin) - 1)); % compute p val
        clear tmat;
    end
    
    % write out results
    system(['echo Writing out results >> ' log_file]);
    disp('Writing out results');


    for i = 1:length(seed_mask_index)

        current_index = seed_mask_index(i);
        [y,x,z] = ind2sub(size(seed_mask.vol), current_index);
        y = y - 1;
        x = x - 1;
        z = z - 1;

        corr_filename = fullfile(output_dir, ['lh.' num2str(x) '_' num2str(y) '_' num2str(z) '.corr.nii.gz']);
        output_mask.vol = reshape(corr_mat(1:lh_num_vertices, i), size(output_mask.vol));
        MRIwrite(output_mask, corr_filename);

        corr_filename = fullfile(output_dir, ['rh.' num2str(x) '_' num2str(y) '_' num2str(z) '.corr.nii.gz']);
        output_mask.vol = reshape(corr_mat(lh_num_vertices+1:end, i), size(output_mask.vol));
        MRIwrite(output_mask, corr_filename);
        
        if(pval > 0)
            prob_filename = fullfile(output_dir, ['lh.' num2str(x) '_' num2str(y) '_' num2str(z) '.mlog10p.nii.gz']);
            output_mask.vol = reshape(pmat(1:lh_num_vertices, i), size(output_mask.vol));
            MRIwrite(output_mask, prob_filename);

            prob_filename = fullfile(output_dir, ['rh.' num2str(x) '_' num2str(y) '_' num2str(z) '.mlog10p.nii.gz']);
            output_mask.vol = reshape(pmat(lh_num_vertices+1:end, i), size(output_mask.vol));
            MRIwrite(output_mask, prob_filename);
        end
    end
end

exit
