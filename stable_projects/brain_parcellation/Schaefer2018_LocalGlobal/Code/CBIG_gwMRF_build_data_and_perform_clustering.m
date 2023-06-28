function [results,prams]=CBIG_gwMRF_build_data_and_perform_clustering(input_fullpaths,output_path,start_idx,end_idx,...
    num_left_cluster,num_right_cluster,smoothcost,num_iterations,num_runs,start_gamma,exponential, iter_reduce_gamma)

    % This function reads in the data, computes the multiplication matrix and performs the clustering
    % Note that reading in the data and building the matrices might take some time and need some disk space and memory
    % If the matrix files have been already created(and still exist) this process will be skiped the next time
    % This function does not give you access to all the parameters of test_graph_clustering. You will need to call
    % test_graph_cut_clustering directly to have all options. See also alex_set_prams.m for parameters
    %
    % Input:
    %   - input_fullpaths = a file containing full paths to all subjects' surf data; 
    %     each line represents a subject with different runs
    %   - output_path = path to where output files are written
    %   - start_idx = for selecting subsets in the subject list, this will be the start index
    %   - Note that start_idx is also used for specifying the starting seed.
    %   - end_idx = for selecting subsets in the subject list, this will be the end index
    %   - num_left_cluster =  number of cluster on the left hemisphere;
    %   - num_righ_cluster = number of cluster for the right hemisphere;
    %   - smoothcost = weight for the smoothness in the MRF, see the paper for more details;
    %   - num_iterations = number of iterations for each random initializations;
    %   - num_runs = number of random initializations
    %
    % Output:
    %   - results = struct containing the clustering results
    %   - prams = struct containing the parameters used in computing the clustering results
    %
    % Example:
    %   - [results,prams]=CBIG_gwMRF_build_data_and_perform_clustering('<your_input_fullpaths_file>',
    %     '<your_output_path>',1,2,50,50,5000,7,2,50000000,15);
    %
    % Written by A. Schaefer and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
    addpath('./lib/');
    
    % hardcode to have less parameters, you may change this to the mesh you used
    fsaverage='fsaverage6';
    
    % create folders
    mkdir(output_path);
    mkdir([output_path,'/time_data/']);
    mkdir([output_path,'/mult_mat/']);
    
    %% actuall computation
    % skip if input data already exists
    if ~(exist([output_path,'/mult_mat/rh_mult_matrix.mat'], 'file') >0)
        [lh_output_file,rh_output_file]=CBIG_gwMRF_build_time_matrix(input_fullpaths,[output_path,'/time_data/'],...
        start_idx,end_idx,fsaverage,'lh_time_matrix.mat','rh_time_matrix.mat');

        [lh_output_mult_mat_file,rh_output_mult_mat_file,dim]=CBIG_gwMRF_build_prod_matrix(lh_output_file,...
        rh_output_file,[output_path,'/mult_mat/'],'lh_mult_matrix.mat','rh_mult_matrix.mat');
    else
        lh_output_mult_mat_file=[output_path,'/mult_mat/lh_mult_matrix.mat'];
        rh_output_mult_mat_file=[output_path,'/mult_mat/rh_mult_matrix.mat'];
        load(rh_output_mult_mat_file,'dim'); % just to get the dimensionality
    end
    
    % this parameter was added for unit test only, therefore it is an optional input argument
    % if not specified, we will assign it with a default value
    if(~exist('iter_reduce_gamma','var'))
        iter_reduce_gamma = 300;
    end

    %% do clustering, there is an internal check if results already exist
    mkdir([output_path,'/clustering/']);
    [prams,results] = CBIG_gwMRF_graph_cut_clustering('lh_avg_file',lh_output_mult_mat_file, ...
    'rh_avg_file',rh_output_mult_mat_file,'left_cluster',num_left_cluster,'right_cluster',num_right_cluster,...
    'iterations',num_iterations,'smoothcost',smoothcost,'start_index',start_idx, 'runs', num_runs,...
    'output_folder' ,[output_path,'/clustering/'],'dim',dim,'exponential',exponential,...
    'start_gamma',start_gamma, 'iter_reduce_gamma', iter_reduce_gamma);

    rmpath('./lib/');
end

