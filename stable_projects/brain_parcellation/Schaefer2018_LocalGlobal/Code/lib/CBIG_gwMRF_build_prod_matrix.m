function [lh_output_file,rh_output_file,dim]=CBIG_gwMRF_build_prod_matrix(lh_input_filename,rh_input_filename,output_path,lh_output_filename,rh_output_filename)

    % This function precomputes the products of timeseries data.
    % The time series data can be created with alex_build_time_matrix
    %
    % Input
    %
    %   - lh_input_filename = left hemisphere matlab timeseries file, created by CBIG_create_time_matrix
    %   - rh_input_filename = right hemisphere matlab timeseries file, created by CBIG_create_time_matrix
    %   - output_path = folder to save the product timeseries matrix
    %   - lh_output_filename = filename for left hemisphere product matrix
    %   - rh_output_filename = filename for right hemisphere product matrix
    % 
    % Output
    %   - lh_output_file = left hemisphere precomputed matrix file
    %   - rh_output_file = right hemisphere precomputed matrix file
    %   - dim = dimension of the matrix file
    % Example
    %   - [lh_output_mult_mat_file,rh_output_mult_mat_file,dim]=CBIG_gwMRF_build_prod_matrix(lh_output_file,rh_output_file,[output_path,'/mult_mat/'],'lh_mult_matrix.mat','rh_mult_matrix.mat');
    %
    % Written by Alexander Schaefer and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
    
    %% left hemisphere
    load(lh_input_filename);
    time_mat=lh_time_mat;
    clear lh_vol;
    for i=1:size(time_mat,1) %% zero mean,make unit norm
        time_mat(i,:)=time_mat(i,:)-mean(time_mat(i,:));
        time_mat(i,:)=time_mat(i,:)/sqrt(sum(time_mat(i,:).^2,2));
    end
    cov_mat=time_mat*time_mat';
    dim=size(time_mat,2);
    clear time_mat
    
    % Save results 
    lh_output_file=[output_path,'/',lh_output_filename];
    save(lh_output_file,'cov_mat','dim','-v7.3');
     
    
    %% right hemisphere  
    load(rh_input_filename);
    time_mat=rh_time_mat;
    clear rh_vol;
    for i=1:size(time_mat,1) %% zero mean,make unit norm
        time_mat(i,:)=time_mat(i,:)-mean(time_mat(i,:));
        time_mat(i,:)=time_mat(i,:)/sqrt(sum(time_mat(i,:).^2,2));
    end
    cov_mat=time_mat*time_mat';
    dim=size(time_mat,2);
    clear time_mat
    
    % Save results 
    rh_output_file=[output_path,'/',rh_output_filename];
    save(rh_output_file,'cov_mat','dim','-v7.3');
     
end
