function CBIG_MMLDA_add_text_to_forest_plot(forest_plot, file_pvalue, out_dir, p_thresh, k, compare_factor_names)
% This function takes in the raw forest plot and p values, and add text to 
% the forest plot, such as 'F2 - F1', 'Latter is worse'.
% The final figure will be saved in the same directory as out_dir, with a
% postfix name of '_compPlot.png' to the original file name.
% Currently only available for k = 2, 3, 4 
%
% Input:
%      - forest_plot:
%        The file name of the raw forest plot
%      - file_pvalue:
%        The file name of the text file containing number of subjects and 
%        p-values
%        file_pvalue is assumed to have the following format: 1st row is the 
%        number of subjects, 2nd row is the omnibus p-value, 3rd to the
%        last rows are the pairwise p-values.
%      - out_dir:
%        The full path to directory where the raw GLM plot & p value text 
%        files are
%      - p_thresh:
%        The final threshold after correcting multiple comparisons
%        p-values stronger than p_thresh will be highlighted in blue.
%      - k:
%        Number of factors
%      - compare_factor_names:
%        (optional) cell array of factor names. Default is {'F2 - F1'}
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% read in the forest plot
img = imread(forest_plot);
[file_path, file_name, ext] = fileparts(forest_plot);

%% draw composite plots
% k = 2
if k == 2
    outerBox = imread('/templates_plots/box_template_k2.png');
    
    img = imresize(img, [80,510]);
    
    mask = outerBox;
    
    for i = 50:(size(img,1)+49)
        for j = 210:635
            if mask(i,j,1)== 255 && mask(i,j,2)==255 && mask(i,j,3)==255
                mask(i,j,:) = img(i-49,j-160,:);
            end
        end
    end
    
    imshow(mask);
    
    % insert N & p values
    fid = fopen(file_pvalue,'r');
    p = fscanf(fid,'%d\n%e\n');
    fclose(fid);
    
    pVals = cell(size(p));
    pVals{1} = num2str(p(1));
    for i = 2:length(p)
        pVals{i} = CBIG_MMLDA_pvalue2str(p(i));
    end
    
    FontOptions = cell(length(pVals),2); % FontWeight & Color
    FontOptions(1,:) = {'normal','black'};
    if p(2) < p_thresh
        FontOptions{i,1} = 'bold';
        FontOptions{i,2} = 'blue';
    else
        FontOptions{i,1} = 'normal';
        FontOptions{i,2} = 'black';
    end
     
    text_pVals = cell(length(pVals),1);
    text_pVals{1} = ['N = ' pVals{1}];
    text_pVals{2} = pVals{2};
    
    position_pVals = [70 70; 670 70]; % Positions to write pVals
    
    for i = 1:length(text_pVals)
        text(position_pVals(i,1),position_pVals(i,2),text_pVals{i},'FontSize',25,'FontWeight',FontOptions{i,1},'Color',FontOptions{i,2},'HorizontalAlignment','center');
    end
    
    % insert labels
    if strfind(forest_plot,'AGE')
        text_labels = {'Latter is older'; 'Former is older'; 'p'; 'Contrast'};
    elseif strfind(forest_plot, 'SEX')
        text_labels = {'Latter is female'; 'Former is female'; 'p'; 'Contrast'};
    elseif strfind(forest_plot, 'IQ')
        text_labels = {'Latter is better'; 'Former is better'; 'p'; 'Contrast'};
    else
        text_labels = {'Latter is worse'; 'Former is worse'; 'p'; 'Contrast'};
    end
    
    position_labels = [225 30; 615 30; 670 30; 430 140]; % Pixel positions to write labels
    AlignOptions = {'left','right','center','center'};
    for i = 1:length(text_labels)
        text(position_labels(i,1),position_labels(i,2),text_labels{i},'FontSize',20,'HorizontalAlignment',AlignOptions{i});
    end
    
    % insert factor names
    if (~exist('compare_factor_names'))
        compare_factor_names = {'F2 - F1'};
    end

    position_fNames = [180 70]; % Pixel positions to write factor names
    for i = 1:length(compare_factor_names)
        text(position_fNames(i,1),position_fNames(i,2),compare_factor_names{i},'FontSize',25,'HorizontalAlignment','center');
    end
    
    % write texts into the image and save it in intputDir
    img_comp = getframe(gcf);
    
    file_name = [file_name '_compPlot.png'];
    imwrite(img_comp.cdata,[out_dir file_name]);
    
    set(gcf,'PaperPositionMode','auto');
    print(gcf, [out_dir file_name],'-dpng','-r0');
    
% k = 3
elseif k == 3
    outerBox = imread('/templates_plots/box_template_k3.png');
    
    img = imresize(img, [240,785]);
    
    mask = outerBox;
    
    for i = 61:(size(img,1)+52)
        for j = 325:960
            if mask(i,j,1)== 255 && mask(i,j,2)==255 && mask(i,j,3)==255
                mask(i,j,:) = img(i-53,j-240,:);
            end
        end
    end
    
    imshow(mask);
    
    
    % insert N & p values
    fid = fopen(file_pvalue,'r');
    p = fscanf(fid,'%d\n%e\n%e\n%e\n%e\n');
    fclose(fid);
    
    pVals = cell(size(p));
    pVals{1} = num2str(p(1));
    for i = 2:length(p)
        pVals{i} = CBIG_MMLDA_pvalue2str(p(i));
    end
    
    FontOptions = cell(length(pVals),2); % FontWeight & Color
    for i = 1:2
        FontOptions(i,:) = {'normal','black'};
    end
    for i = 3:length(pVals)
        if p(i) < p_thresh
            FontOptions{i,1} = 'bold';
            FontOptions{i,2} = 'blue';
        else
            FontOptions{i,1} = 'normal';
            FontOptions{i,2} = 'black';
        end
    end
    
    text_pVals = cell(length(pVals),1);
    text_pVals{1} = ['N = ' pVals{1}];
    text_pVals{2} = ['p = ' pVals{2}];
    
    for i = 3:length(pVals)
        text_pVals{i} = pVals{i};
    end
    
    position_pVals = [110 130; 110 190; 1030 100; 1030 170; 1030 230]; % Positions to write pVals
    
    for i = 1:length(text_pVals)
        text(position_pVals(i,1),position_pVals(i,2),text_pVals{i},'FontSize',28,'FontWeight',FontOptions{i,1},'Color',FontOptions{i,2},'HorizontalAlignment','center');
    end
    
    % insert labels
    if strfind(forest_plot,'AGE')
        text_labels = {'Latter is older'; 'Former is older'; 'p'; 'Contrast'};
    elseif strfind(forest_plot, 'SEX')
        text_labels = {'Latter is female'; 'Former is female'; 'p'; 'Contrast'};
    elseif strfind(forest_plot, 'IQ')
        text_labels = {'Latter is better'; 'Former is better'; 'p'; 'Contrast'};
    elseif strfind(forest_plot, 'ICV')
        text_labels = {'Latter is more atrophied'; 'Former is more atrophied'; 'p'; 'Contrast'};
    else
        text_labels = {'Latter is worse'; 'Former is worse'; 'p'; 'Contrast'};
    end
    
    position_labels = [345 45; 950 45; 1030 45; 650 330]; % Pixel positions to write labels
    AlignOptions = {'left','right','center','center'};
    for i = 1:length(text_labels)
        text(position_labels(i,1),position_labels(i,2),text_labels{i},'FontSize',22,'HorizontalAlignment',AlignOptions{i});
    end
    
    % insert factor names
    if (~exist('compare_factor_names'))
        compare_factor_names = {'F2 - F1'; 'F3 - F1'; 'F3 - F2'};
    end
    %compare_factor_names = {'L - T';'C - T';'C - L'};
    position_fNames = [270 105; 270 170; 270 230]; % Pixel positions to write factor names
    for i = 1:length(compare_factor_names)
        text(position_fNames(i,1),position_fNames(i,2),compare_factor_names{i},'FontSize',28,'HorizontalAlignment','center');
    end
    
    % write texts into the image and save it in intputDir
    img_comp = getframe(gcf);
    
    file_name = [file_name '_compPlot.png'];
    imwrite(img_comp.cdata,[out_dir file_name]);
    
    set(gcf,'PaperPositionMode','auto');
    print(gcf, [out_dir file_name],'-dpng','-r0');
 
% k = 4
else
    outerBox = imread('/templates_plots/box_template_k4.png');
    
    img = imresize(img, [300,520]);
    
    mask = outerBox;
    
    for i = 61:(size(img,1)+39)
        for j = 160:(size(img,2)+154)
            if mask(i,j,1)== 255 && mask(i,j,2)==255 && mask(i,j,3)==255
                mask(i,j,:) = img(i-40,j-155,:);
            end
        end
    end
    
    imshow(mask);

    % insert N & p values
    fid = fopen(file_pvalue,'r');
    p = fscanf(fid,'%d\n%e\n%e\n%e\n%e\n%e\n%e\n%e\n');
    fclose(fid);
    
    pVals = cell(size(p));
    pVals{1} = num2str(p(1));
    for i = 2:length(p)
        pVals{i} = CBIG_MMLDA_pvalue2str(p(i));
    end
    
    FontOptions = cell(length(pVals),2); % FontWeight & Color
    for i = 1:2
        FontOptions(i,:) = {'normal','black'};
    end
    for i = 3:length(pVals)
        if p(i) < p_thresh
            FontOptions{i,1} = 'bold';
            FontOptions{i,2} = 'blue';
        else
            FontOptions{i,1} = 'normal';
            FontOptions{i,2} = 'black';
        end
    end
    
    text_pVals = cell(length(pVals),1);
    text_pVals{1} = ['N = ' pVals{1}];
    text_pVals{2} = ['p = ' pVals{2}];
    
    for i = 3:length(pVals)
        text_pVals{i} = pVals{i};
    end
    
    position_pVals = [70 170; 70 200; 670 85; 670 125; 670 165; 670 205; 670 245; 670 285]; % Positions to write pVals
    
    for i = 1:length(text_pVals)
        text(position_pVals(i,1),position_pVals(i,2),text_pVals{i},'FontSize',24,'FontWeight',FontOptions{i,1},'Color',FontOptions{i,2},'HorizontalAlignment','center');
    end
    
    % insert labels
    if strfind(forest_plot,'AGE')
        text_labels = {'Latter is older'; 'Former is older'; 'p'; 'Contrast'};
    elseif strfind(forest_plot, 'SEX')
        text_labels = {'Latter is female'; 'Former is female'; 'p'; 'Contrast'};
    elseif strfind(forest_plot, 'IQ')
        text_labels = {'Latter is better'; 'Former is better'; 'p'; 'Contrast'};
    elseif strfind(forest_plot, 'ICV')
        text_labels = {'Latter is more atrophied'; 'Former is more atrophied'; 'p'; 'Contrast'};
    else
        text_labels = {'Latter is worse'; 'Former is worse'; 'p'; 'Contrast'};
    end
    
    position_labels = [230 50; 620 50; 670 50; 420 350]; % Pixel positions to write labels
    AlignOptions = {'left','right','center','center'};
    for i = 1:length(text_labels)
        text(position_labels(i,1),position_labels(i,2),text_labels{i},'FontSize',20,'HorizontalAlignment',AlignOptions{i});
    end
    
    % insert factor names
    if (~exist('compare_factor_names'))
        compare_factor_names = {'F2 - F1'; 'F3 - F1'; 'F4 - F1'; 'F3 - F2'; 'F4 - F2'; 'F4 - F3'};
    end

    position_fNames = [180 85; 180 125; 180 165; 180 205; 180 245; 180 285]; % Pixel positions to write factor names
    for i = 1:length(compare_factor_names)
        text(position_fNames(i,1),position_fNames(i,2),compare_factor_names{i},'FontSize',25,'HorizontalAlignment','center');
    end
    
    % write texts into the image and save it in intputDir
    img_comp = getframe(gcf);
    
    file_name = [file_name '_compPlot.png'];
    imwrite(img_comp.cdata,[out_dir file_name]);
    
    set(gcf,'PaperPositionMode','auto');
    print(gcf, [out_dir file_name],'-dpng','-r0');
    
end

end



