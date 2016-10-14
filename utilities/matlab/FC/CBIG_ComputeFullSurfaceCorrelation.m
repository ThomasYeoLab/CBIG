function CBIG_ComputeFullSurfaceCorrelation(output_file, varargin_text1, varargin_text2, pval)

% CBIG_ComputeFullSurfaceCorrelation(output_file, varargin_text1, varargin_text2, pval)
% This function is used to compute full surface correlation
% [lh_hemi rh_hemi] by [lh_hemi rh_hemi]. 
% Input:
%      -output_file: output .mat file.
%      -varargin_text1: left hemisphere data list, each row corresponds to
%       a .nii/.nii.gz fMRI data file.
%      -varargin_text2: right hemisphere data list, each row corresponds to
%       a .nii/.nii.gz fMRI data file.
%      -pval: If pval > 0 then this function will compute statistics
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% read in both left and right text files.
fid = fopen(varargin_text1, 'r');
i = 0;
while(1);
    tmp = fgetl(fid);
    if(tmp == -1)
        break
    else
        i = i + 1;
        varargin1{i} = tmp;
    end
end
fclose(fid);

% read in both left and right text files.
fid = fopen(varargin_text2, 'r');
i = 0;
while(1);
    tmp = fgetl(fid);
    if(tmp == -1)
        break
    else
        i = i + 1;
        varargin2{i} = tmp;
    end
end
fclose(fid);

if(nargin < 4)
    pval = 0;
else
    if(ischar(pval))
        pval = str2num(pval);
    end
end

log_file = [output_file '.log'];
delete(log_file);

% Compute correlation
for i = 1:length(varargin1)

    disp(num2str(i));
    system(['echo ' num2str(i) ' >> ' log_file]);

    C1 = textscan(varargin1{i}, '%s');
    C1 = C1{1};

    C2 = textscan(varargin2{i}, '%s');
    C2 = C2{1};

    for j = 1:length(C1)
        input = C1{j};
        input_series = MRIread(input);
        t_series1 = single(transpose(reshape(input_series.vol, size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), size(input_series.vol, 4))));

        input = C2{j};
        input_series = MRIread(input);
        t_series2 = single(transpose(reshape(input_series.vol, size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), size(input_series.vol, 4))));

        t_series = [t_series1 t_series2];

        % normalize series (note that series are now of dimensions: T x N)
        t_series = bsxfun(@minus, t_series, mean(t_series, 1));
        t_series = bsxfun(@times, t_series, 1./sqrt(sum(t_series.^2, 1)));

        corr_mat = t_series' * t_series;
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
output = output / length(varargin1);
corr_mat = tanh(output);
disp(['isnan: ' num2str(sum(isnan(corr_mat(:)))) ' out of ' num2str(numel(corr_mat))]);
corr_mat(isnan(corr_mat)) = 0;
clear output;

if(pval > 0)
    system(['echo Computing statistics >> ' log_file]);
    disp('Computing statistics');

    zmat = zmat/length(varargin1); % compute mean
    zsqmat = zsqmat/length(varargin1);
    stdz = zsqmat - zmat.^2; stdz(stdz < 0) = 0;
    stdz = sqrt(stdz)*length(varargin1)/(length(varargin1) - 1); %compute std
    clear zsqmat;
    tmat = sqrt(length(varargin1))*zmat./(stdz + eps); % compute t stats
    clear stdz; clear zmat;
    pmat = -log10(tcdf(-abs(tmat), length(varargin1) - 1)); % compute p val
    clear tmat;
end

% write out results
system(['echo Writing out results >> ' log_file]);
disp('Writing out results');

if(pval > 0) 
    save(output_file, 'corr_mat', 'pmat', '-v7.3');
else
    save(output_file, 'corr_mat', '-v7.3');
end
exit
