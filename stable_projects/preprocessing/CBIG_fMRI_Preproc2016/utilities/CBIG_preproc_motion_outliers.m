function CBIG_preproc_motion_outliers(DVARS_file,FD_metric_file,FD_th,DV_th,discard_seg,FD_metric,output_dir_name)

% CBIG_preproc_motion_outliers(DVARS_file,FD_metric_file,FD_th,DV_th,discard_seg,FD_metric,output_dir_name)
%
% This function is used to find the outliers with DVARS threshold DV_th and
% FDRMS/FD threshold FD_th. After thresholding, one frame before two frame after
% will also be removed. The default DV_th is 50, default FD_th is 0.2
% INPUT:
%   - DVARS_file:
%       a file contains DVARS values as a TX1 column, T is the number of
%       frames (include the first frame, whose value supposed to be 0)
%
%   - FD_metric_file:
%       a file contains FDRMS/FD values as a TX1 column, T is the number of
%       frames (include the first frame, whose value supposed to be 0)
%
%   - FD_th:
%       the threshold of FDRMS/FD, default is 0.2
%
%   - DV_th:
%       the threshold of DVARS, default is 50
%
%   - discard_seg:
%       uncensored segments of data lasting fewer than discard_seg
%       contiguous frames will be removed, default is 5.
%
%   -FD_metric:
%       A string to determine the FD metric used.
%       Can contain only either 'FD' or 'FDRMS'.
%
%   - output_dir_name:
%       the directory to store all output figures (full path).
%
% OUTPUT:
%   txt file: SUBJECT_bldXXX_<FD_metric>xxx_DVARSxxx_motion_outliers.txt, this file contains a Tx1
%             vector where 1 means kept, 0 means removed.
%   Figure: SUBJECT_bldXXX_<FD_metric>xxx_DVARSxxx.png
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



%% Read in file and construct DVARS. FDRMS/FD vectors
fd = dlmread(FD_metric_file);
dvars = dlmread(DVARS_file);

%% Check input variables
if size(fd,2) > 1
    error('Input argument ''FD_metric_file'' should list values in columns');
end

if size(dvars,2) > 1
    error('Input argument ''DVARS_file'' should list values in columns');
end

if ~strcmp(FD_metric, 'FD') && ~strcmp(FD_metric, 'FDRMS')
    error('ERROR: FD_metric can only be FDRMS or FD.');
end

%% Detect outliers with DVARS threshold
DV_th = str2num(DV_th);
% 1 means the frame below the threshold is kept, 0 means the frame above the threshold is removed
DV_vect = (dvars <= DV_th);

%% Detect outliers with FDRMS/FD threshold
FD_th = str2num(FD_th);
% 1 means the frame below the threshold is kept, 0 means the frame above the threshold is removed
FD_vect = (fd <= FD_th);

%% Continue removing one frame before and two frames after
DV_censor = (DV_vect & [DV_vect(2:end,1); 1] & [1; DV_vect(1:end-1)] & [1;1; DV_vect(1:end-2)]);
FD_censor = (FD_vect & [FD_vect(2:end,1); 1] & [1; FD_vect(1:end-1)] & [1;1; FD_vect(1:end-2)]);

%% Plot DVARS and FDRMS/FD with threshold line
% Plot FDRMS/FD
pos_set = [0.08 0.55 0.89 0.4; ...
    0.08 0.09 0.89 0.37];

figure;
set(gcf,'Visible','off');
set(gcf, 'Position', [0,0,960,800]);
g = subplot(2,1,1);
set(g, 'Position', pos_set(1,:));
hold on
plot(fd, 'LineWidth', 2);
ti = title([FD_metric ' (in mm)']);
set(ti, 'FontSize', 18);
ylabel('FD (mm)', 'FontSize', 18);
% Get automatically generated x-axis limits
xlim = get(gca, 'XLim');
line([xlim(1) xlim(2)], [FD_th FD_th], 'color', 'r', 'LineWidth', 1.5);
% Let MATLAB automatically determine major and minor ticks
set(gca, 'XMinorTick', 'on');
set(gca, 'FontSize', 18)
hold off
set(gcf, 'PaperPositionMode', 'auto')
clear xlim

% Plot DVARS
g = subplot(2,1,2);
set(g, 'Position', pos_set(2,:));
hold on
plot(dvars, 'LineWidth', 2);
ti = title('DVARS');
set(ti, 'FontSize', 18);
xlabel('Frame Number', 'FontSize', 18);
ylabel('DVARS', 'FontSize', 18);
xlim = get(gca, 'XLim');
line([xlim(1) xlim(2)], [DV_th DV_th], 'color', 'r', 'LineWidth', 1.5);
set(gca, 'XMinorTick', 'on');
set(gca, 'FontSize', 18)
hold off
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, [output_dir_name '_' FD_metric num2str(FD_th) ...
    '_DVARS' num2str(DV_th) '.png'], '-dpng');
clear xlim
close(gcf);

%% Combine DVARS and FDRMS/FD so that a frame is either above FD_th or DV_th, it will be removed
DV_FD_censor = [DV_censor & FD_censor];

%% remove uncensored segments of data lasting fewer than discard_seg contiguous frames
discard_seg = str2num(discard_seg);
DV_FD_censor_segment = diff([0; DV_FD_censor; 0]);
seg_start = find(DV_FD_censor_segment == 1);
seg_end = find(DV_FD_censor_segment == -1)-1;

for idx = 1:length(seg_start)
    if(seg_end(idx)-seg_start(idx) + 1 < discard_seg)
        DV_FD_censor(seg_start(idx):seg_end(idx)) = 0;
    end
end

%% Save out the outlier vector into a text file
outfilename_index = fopen([output_dir_name '_' FD_metric num2str(FD_th) ...
    '_DVARS' num2str(DV_th) '_motion_outliers.txt'], 'w+');
fprintf(outfilename_index, '%d\n', DV_FD_censor);
fclose(outfilename_index);
