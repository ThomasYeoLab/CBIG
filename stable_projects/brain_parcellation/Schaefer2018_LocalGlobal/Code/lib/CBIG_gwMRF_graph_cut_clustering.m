function [prams,results]=CBIG_gwMRF_graph_cut_clustering(varargin)

    % This function performs MRF clustering using Graphcuts; modelling the data as VMF distributions.
    % Results will be stored in the results struct. The used parameters will be saved in the prams struct.

    % Input
    %
    % Possible parameters are (please also check alex_set_prams.m):
    %   - start_index = will be the start seed e.g. 1 or 5982
    %   - runs = will be the number of runs after the seed e.g. 10, together with startseed 1 this would run seed 1, 2, 3, ..., 10
    %   - cluster = will be the number of cluster e.g. 160
    %   - datacost and smoothcost = scalars which give the scaling of the respective costs e.g. 1000 or 0.1
    %   - expo = exponetial k we use on the prior, the formular is exp(-k * gradient) - exp(-k)
    %   - grad_prior = can be 'smoothwm' or 'midthickness'
    %   - output_folder = can be set as a string, will be by default './output/'
    %   - lh_avg_file = is the left hemisphere input profile, can be set as a string, will be by default './input/lh.avg_Train_GSP_release_weighted_fs6_ROIfs4_profiles.nii.gz';
    %   - rh_avg_file = is the right hemisphere input, can be set as a string, will be by default './input/rh.avg_Train_GSP_release_weighted_fs6_ROIfs4_profiles.nii.gz';
    %   - fsaverage = is a string and will be by default 'fsaverage6'
    %   - iterations = is an integer number e.g. 100
    % If parameters are not set they will be set to some default value which can be found in CBIG_gwMRF_set_prams.m
    % Paramters should be given in a form like 'expo','25' or 'expo',25 the order does not matter but the value(25) has to follow the keyword (expo)
    %
    % Output:
    %   - results = struct containing the clustering results
    %   - prams = struct containing the parameters used in computing the clustering results
    % 
    %Example
    %   - [prams,results]=CBIG_gwMRF_graph_cut_clustering('lh_avg_file',lh_output_mult_mat_file,'rh_avg_file',rh_output_mult_mat_file,'left_cluster','50','right_cluster','50','iterations',100,'smoothcost','100000','start_index','1', 'runs', num_runs,'output_folder' ,[output_path,'/clustering/'],'dim',size(time_mat,2));
    %
    %Written by Alexander Schaefer and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    % compute cluster
    prams=CBIG_gwMRF_set_prams(varargin{:});
    fprintf(prams.fileID,'Will iterate from initialization %i to initialization %i \n',prams.start_index,prams.start_index+prams.runs-1)
    for i=prams.start_index:(prams.start_index+prams.runs-1)
        prams.seed=i;
        if(exist([prams.output_folder,prams.output_name,'_seed_',num2str(prams.seed),'_Energy.mat'],'file'))% If file for this seed already exists we just skip
            fprintf(prams.fileID,'File for seed %i already exists, skipping! \n',prams.seed)
            continue
        else
            [results]=CBIG_gwMRF_graph_cut_clustering_iter_split(prams);
            save([prams.output_folder,prams.output_name,'_seed_',num2str(prams.seed),'.mat'],'results','prams')
            Energy=results.E;
            save([prams.output_folder,prams.output_name,'_seed_',num2str(prams.seed),'_Energy.mat'],'Energy') % save energy in a seperate file
        end
    end
    fclose(prams.fileID);%%% close output file
end
