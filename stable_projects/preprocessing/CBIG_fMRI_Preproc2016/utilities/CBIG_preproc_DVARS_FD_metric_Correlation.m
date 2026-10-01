function CBIG_preproc_DVARS_FD_metric_Correlation(DVARS_file,FD_metric_file,FD_metric,output_dir_name)

% CBIG_preproc_DVARS_FD_metric_Correlation(DVARS_file,FD_metric_file,FD_metric,output_dir_name)
%
% This function is used to compute the correlation between DVARS and FDRMS/FD.
% This correlation shows whether DVARS and FDRMS/FD are able to caputure
% the same movement trend. High correlation indicates both DVARS and FDRMS/FD
% caputure the same motion feature.
%
% INPUT:
%      -DVARS_file:
%       a file contains DVARS values as a T x 1 column, where T is the
%       number of frames (including the first frame, whose value is 0)
%
%      -FD_metric_file:
%       a file contains FDRMS/FD values as a T x 1 column, where T is the
%       number of frames (including the first frame, whose value is 0)
%
%      -FD_metric:
%       A string to determine the FD metric used.
%       Can contain only either 'FD' or 'FDRMS'.
%
%      -output_dir_name:
%       A string containing the path to the output directory.
%
% OUTPUT:
%      Figure: SUBJECT_bldXXX_DVARS_<FD_metric>_correlation.png
%
% EXAMPLE:
%     DVARS_file = 'XXX/bold/mc/sub005_bld001_rest_skip8_stc_motion_outliers_DVARS'
%     FD_metric_file = 'XXX/bold/mc/sub005_bld001_rest_skip8_stc_motion_outliers_FDRMS'
%     FD_metric = 'FDRMS'
%     output_dir_name = 'XXX/qc/sub005_bld001'
%     CBIG_preproc_DVARS_FD_metric_Correlation(DVARS_file, FD_metric_file, ...
%         FD_metric, output_dir_name)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Read in file and construct DVARS. FDRMS/FD vectors
fd = csvread(FD_metric_file);
dvars = csvread(DVARS_file);

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

%% Compute the correlation, discard the first frame
coe = CBIG_corr(fd(2:end), dvars(2:end));

%% Plot and save the figure
figure;
set(gcf, 'Visible', 'off');
set(gcf, 'Position', [0,0,960,800]);
ti = title(['Correlation between DVARS and ' FD_metric ': ' num2str(coe)]);
set(ti, 'FontSize', 22);
xlabel(FD_metric, 'FontSize', 22);
ylabel('DVARS', 'FontSize', 22);
set(gca, 'FontSize', 22);
outerpos = get(gca, 'OuterPosition');
tight = get(gca, 'TightInset');
left = outerpos(1) + tight(1);
bottom = outerpos(2) + tight(2);
width = outerpos(3) - tight(1) - tight(3);
height = outerpos(4) - tight(2) - tight(4);
set(gca, 'Position', [left, bottom, width, height])

hold on;
scatter(fd(2:end), dvars(2:end));
%Create the fitting line
p=polyfit(fd(2:end), dvars(2:end),1);
r=polyval(p, fd(2:end));
plot(fd(2:end), r, '-r', 'LineWidth', 2);
hold off;
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, [output_dir_name '_DVARS_' FD_metric '_correlation'], '-dpng');
close(gcf);

