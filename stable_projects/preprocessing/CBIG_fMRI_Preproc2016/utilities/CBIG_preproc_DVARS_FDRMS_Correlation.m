function CBIG_preproc_DVARS_FDRMS_Correlation(DVARS_file,FDRMS_file,output_dir_name)

%CBIG_preproc_DVARS_FDRMS_Correlation(DVARS_file,FDRMS_file,output_dir_name)
%
% This function is used to compute the correlation between DVARS and FDRMS. 
% This correlation shows whether DVARS and FDRMS are able to caputure      
% the same movement trend. High correlation indicates both DVARS and FDRMS 
% caputure the same motion feature.                                            
%
% INPUT:                                                                   
%      -DVARS_file: 
%       a file contains DVARS values as a T x 1 column, where T is the 
%	number of frames (including the first frame, whose value is 0)
%                    
%      -FDRMS_file: 
%       a file contains FDRMS values as a T x 1 column, where T is the 
%	number of frames (including the first frame, whose value is 0)                    
%
% OUTPUT:                                                                 
%      Figure: SUBJECT_bldXXX_DVARS_FDRMS_correlation.jpg                    
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%% Read in file and construct DVARS. FDRMS vectors
fd = csvread(FDRMS_file);
dvars = csvread(DVARS_file);

%% Compute the correlation, discard the first frame
coe = CBIG_corr(fd(2:end), dvars(2:end));

%% Plot and save the figure
figure;
set(gcf, 'Visible', 'off');
set(gcf, 'Position', [0,0,960,800]);
ti = title(['Correlation between DVARS and FDRMS:' num2str(coe)]);
set(ti, 'FontSize', 22);
xlabel('FDRMS', 'FontSize', 22);
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
print(gcf, [output_dir_name '_DVARS_FDRMS_correlation'], '-dpng');
close(gcf);

