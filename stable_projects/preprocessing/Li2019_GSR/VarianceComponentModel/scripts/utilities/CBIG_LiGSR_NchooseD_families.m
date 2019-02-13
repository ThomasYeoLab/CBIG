function CBIG_LiGSR_NchooseD_families( fullsub_list, family_csv, family_header, subject_header, ...
    d, num_samples, outdir, output_prefix )

% CBIG_LiGSR_NchooseD_families( fullsub_list, family_csv, family_header, ...
%     d, num_samples, outdir, output_prefix )
% 
% This function creates a list of subjects to be removed for each delete-d
% jackknife sample.
% 
% Inputs:
%   - fullsub_list
%     A string. The full subject list (full path). Each line is one subject ID.
%    
%   - family_csv
%     A string. The full path of the CSV file that contains the family
%     information of all subjects. When all subjects are independent, you
%     can input 'NONE'.
% 
%   - family_header
%     A string corresponding to the header of the family ID column in
%     "family_csv" file. When all subjects are independent, you can input
%     'NONE'.
% 
%   - subject_header
%     A string that corresponds to the header of the subject ID column in
%     "family_csv" file. When all subjects are independent, you can input
%     'NONE'.
%   
%   - d
%     A scalar or a string, the number of subjects to be removed for each
%     jackknife sample, e.g. 431.
% 
%   - num_samples
%     A scalar or a string, the total number of jackknife samples, e.g.
%     1000.
% 
%   - outdir
%     A string. The full path of the output directory.
%   
%   - output_prefix
%     A string, the prefix to the filename of each output subject list. The
%     output filenames will be
%     [outdir '/' output_prefix '_choose' num2str(d) '_set' num2str(i) '.txt']
%     where i ranges from  1 to num_samples.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

[sub_names, num_sub] = CBIG_text2cell(fullsub_list);

if(ischar(d))
    d = str2num(d);
end

if(ischar(num_samples))
    num_samples = str2num(num_samples);
end

seed = 100;

if(~exist(outdir, 'dir'))
    mkdir(outdir);
end

% check if all lists of jackknife samples exist
exist_flag = 1;
for i = 1:num_samples
    if (~exist(fullfile(outdir, [output_prefix '_choose' num2str(d) '_set' num2str(i) '.txt']), 'file'))
        exist_flag = 0;
        break
    end
end

if(exist_flag==0)
    rng(seed);
    for i = 1:num_samples
        if(~strcmpi(family_csv, 'none') && isempty(family_csv))
            % if subjects are related
            family_IDs = CBIG_parse_delimited_txtfile(family_csv, {family_header}, [], ...
                subject_header, subjects, delimiter);
            unique_families = unique(family_IDs);
            num_fam = length(unique_families);
            fprintf('Total number of families: %d.\n', num_fam);
            
            sub_perfam = cell(num_fam, 1);
            for f = 1:num_fam
                fam_ind = strcmp(family_IDs, unique_families{f})==1;
                sub_perfam{f} = sub_names(fam_ind)';
            end
            fam_chosen = datasample(sub_perfam, d, 2, 'Replace', false);
            fam_chosen = fam_chosen(:);
            fam_chosen(cellfun(@isempty, fam_chosen)) = [];
            CBIG_cell2text(fam_chosen, fullfile(outdir, ...
                [output_prefix '_choose' num2str(d) '_set' num2str(i) '.txt']));
        else
            % if subjects are unrelated
            sub_chosen = datasample(sub_names, d, 'Replace', false);
            CBIG_cell2text(sub_chosen, fullfile(outdir, ...
                [output_prefix '_choose' num2str(d) '_set' num2str(i) '.txt']));
        end
    end
else
    warning('All subject lists exist. Skipping ...\n')
end

end

