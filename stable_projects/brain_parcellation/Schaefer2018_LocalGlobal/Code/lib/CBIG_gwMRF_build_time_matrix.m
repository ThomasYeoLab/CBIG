function [lh_output_file,rh_output_file]=CBIG_gwMRF_build_time_matrix(input_fullpaths,output_path,start_idx,end_idx,fsaverage,lh_output_file,rh_output_file)

    % This function concatenates timeseries data from several subjects
    % Possible inputs for our current setting:
    %
    % Input
    %
    %   - input_fullpaths = a file containing full paths to all subjects' surf data; each line represents a subject with different runs
    %   - output_path = path to where output files are written
    %   - start_idx = for selecting subsets in the subject list, this will be the start index
    %   - end_idx = for selecting subsets in the subject list, this will be the end index
    %   - fsaverage = which fsaverage resolution to use
    %   - lh_output_filename = filename for left hemisphere product matrix
    %   - rh_output_filename = filename for right hemisphere product matrix
    %
    % Ouput
    %   - lh_output_file = matrix containing the concatenated timeseries data of the left hemisphere
    %   - rh_output_file = matrix containing the concatenated timeseries data of the right hemisphere
    %
    % Example
    %   - [lh_output_file,rh_output_file]=CBIG_gwMRF_build_time_matrix(input_fullpaths,[output_path,'/time_data/'],start_idx,end_idx,'fsaverage6','lh_time_matrix.mat','rh_time_matrix.mat');
    %   - Example for data used in the paper
    %       - Release Data
    %           - input_path='$HOME/projects/MRF/data/GSP_official/subject_path.txt'; 
    %       - Test
    %           - input_filename='$HOME/projects/MRF/data/GSP_official/GSP_CSV/test/test_filenames_corrected.csv'
    %           - input_subjectname='$HOME/projects/MRF/data/GSP_official/GSP_CSV/test/test_subjectname_corrected.csv'
    %       - Train
    %           - input_filename='$HOME/projects/MRF/data/GSP_official/GSP_CSV/train/train_filenames_corrected.csv'
    %           - input_subjectname='$HOME/projects/MRF/data/GSP_official/GSP_CSV/train/train_subjectnames_corrected.csv'
    %       - Full
    %           - input_filename='$HOME/projects/MRF/data/GSP_official/GSP_CSV/full/full_filenames_corrected.csv'
    %           - input_subjectname='$HOME/projects/MRF/data/GSP_official/GSP_CSV/full/full_subjectnames_corrected.csv'
    %
    % Written by Alexander Schaefer and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    if(ischar(start_idx)) 
        start_idx=str2num(start_idx);
    end

    if(ischar(end_idx)) 
        end_idx=str2num(end_idx);
    end

    full_paths=table2cell(readtable(input_fullpaths,'Delimiter',' ','ReadVariableNames',false));
    [num_subs, num_scans]=size(full_paths);
    
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', fsaverage, 'inflated', 'cortex'); %for initializing group average, which fsaverage should be given from outside
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', fsaverage, 'inflated', 'cortex');
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%to preallocate matrix size
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    matrix_number_of_scans=0;% to preallocate matrix size
    
    for k=start_idx:end_idx
        for i=1:num_scans
                if (~isnan(full_paths{k,i}))
                    if(~isempty(full_paths{k,i}))
                        matrix_number_of_scans=matrix_number_of_scans+1;
                    end
                end
        end
    end
    addr=full_paths{1,1};
    hemi=MRIread(addr);
    length_of_scan=size(hemi.vol,4);
    lh_time_mat=zeros([size(find(lh_avg_mesh.MARS_label==2),2),length_of_scan*matrix_number_of_scans],'single');%initialize group average
    rh_time_mat=zeros([size(find(rh_avg_mesh.MARS_label==2),2),length_of_scan*matrix_number_of_scans],'single');%initialize group average
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_subjects=0;
    matrix_number_of_scans=0;
    for k=start_idx:end_idx
        fprintf('It is subject number %g ',k);
        sub_num_scans=0;
        
        for l=1:num_scans % determine subject number of scans
            if(length(full_paths{k,l})>4) % assume folder and filename are longer than 4 chars
                sub_num_scans=sub_num_scans+1;
            end
        end
        num_subjects=num_subjects+1;
        scans(k)=sub_num_scans;
        
        for i=1:num_scans
                if (~isnan(full_paths{k,i}))
                    if(~isempty(full_paths{k,i}))
                        addr=full_paths{k,i};
                        lh_input=addr;
                        hemi_index=strfind(lh_input, 'lh');
                        rh_input=lh_input;
                        rh_input(hemi_index:hemi_index+1)='rh';
                        files_used_lh{i}=lh_input;%%just as a log
                        files_used_rh{i}=rh_input;
                        matrix_number_of_scans=matrix_number_of_scans+1;
                        for j=1:2% loop through hemispheres
                            if (j==1)
                                hemi=MRIread(lh_input);
                                vol=reshape(hemi.vol,[size(hemi.vol,1)*size(hemi.vol,2)*size(hemi.vol,3) size(hemi.vol,4)]);
                                vol=single(vol(lh_avg_mesh.MARS_label==2,:));
                                lh_input
                            elseif(j==2)
                                hemi=MRIread(rh_input);
                                vol=reshape(hemi.vol,[size(hemi.vol,1)*size(hemi.vol,2)*size(hemi.vol,3) size(hemi.vol,4)]);
                                vol=single(vol(rh_avg_mesh.MARS_label==2,:));                               
                                rh_input
                            end
                                vol=bsxfun(@minus,vol,mean(vol,2)); % row wise demeaning
                                vol=bsxfun(@rdivide,vol,std(vol',1)'); % row wise standardization
                                
                            if (j==1)
                                lh_time_mat(:,(matrix_number_of_scans-1)*length_of_scan+1:(matrix_number_of_scans)*length_of_scan)=vol;
                            elseif(j==2)
                                rh_time_mat(:,(matrix_number_of_scans-1)*length_of_scan+1:(matrix_number_of_scans)*length_of_scan)=vol;
                            end
                                
                        end
                    end
                end
        end
    end
     %% Save results 
     lh_output_file=[output_path,'/',lh_output_file,];
     save(lh_output_file,'lh_time_mat','scans','files_used_lh','-v7.3');
     rh_output_file=[output_path,'/',rh_output_file];
     save(rh_output_file,'rh_time_mat','scans','files_used_rh','-v7.3');
end
