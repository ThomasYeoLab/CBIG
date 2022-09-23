function CBIG_DiffProc_mrtrix_fill_na(log_file)
% This function reads the MRtrix log file. Entries in the connectome that correspond to nodes that are reported 
% to be missing from the parcellation image are replaced with NaN from 0. The function reads the path to the 
% connectomes directly from the log file.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    fid = fopen(log_file, 'r');
    line = fgetl(fid);
    while ischar(line)
        % get current parcellation
        if contains(line,'Current connectome')
            line_split = strsplit(line);
            curr_parc = line_split{6};
            SC_matrix = csvread(curr_parc);
        elseif strcmp(line,'tck2connectome: [WARNING] The following nodes are missing from the parcellation image:')   
            line = fgetl(fid);
            % get missing nodes
            line_split = strsplit(line, '] ');
            nodes = strsplit(line_split{2},', ');
            nodes = cellfun(@str2num,nodes);
            % replace with NaN
            SC_matrix(:,nodes) = NaN;
            SC_matrix(nodes,:) = NaN;
            csvwrite(curr_parc,SC_matrix)
        end
        line = fgetl(fid);
    end
    fclose(fid);

end
