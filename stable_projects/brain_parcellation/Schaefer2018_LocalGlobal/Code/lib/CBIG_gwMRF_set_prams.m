function prams=CBIG_gwMRF_set_prams(varargin)

    % This function sets the prams file, either from default values or from the input varargin
    % some variables and their explenation and/or possible setting. 
    % If not other specified 1 means true and 0 means False.
    %
    % 
    % Input
    %
    %   - start_index - will be the start seed e.g. 1 or 5982
    %   - runs - the number of runs after the seed e.g. 10, 
    %   together with startseed 1 this would run seed 1, 2, 3, ..., 10
    %   - cluster - the number of cluster e.g. 160
    %   - datacost and smoothcost - scalars which give the scaling of the respective costs e.g. 1000 or 0.1
    %   - expo - exponetial k we use on the prior, the formular is exp(-k * gradient) - exp(-k)
    %   - grad_prior - this is gradient prior for the parcellation. the paper
    %   - used option 'gordon'. You could define your own gradient priors
    %   - output_folder - can be set as a string, will be by default './output/'
    %   - lh_avg_file - is the left hemisphere input profile, can be set as a string, will be by default 
    %   './lib/input/lh.avg_Train_GSP_release_weighted_fs6_ROIfs4_profiles.nii.gz'; 
    %   - rh_avg_file - is the right hemisphere input, can be set as a string, will be by default 
    %   './lib/input/rh.avg_Train_GSP_release_weighted_fs6_ROIfs4_profiles.nii.gz';
    %   - fsaverage - is a string and will be by default 'fsaverage6'
    %   - iterations - is an integer number e.g. 100
    %
    % Output
    %
    %   - prams = a struct containing all the saved parameters
    %
    % Example
    %   - Paramters should be given in a form like 'expo','25' or 'expo',25 
    %   the order does not matter but the value(25) has to follow the keyword (expo)
    %
    % Written by A. Schaefer and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
    p=inputParser;

    %%default values
    defaultseed=5392; 
    % some random number, this seed number is relevant for the random number generator that affects the initalizations
    defaultk=20; % if you dont have different cluster numbers for left and right, this option should work.
    default_left_k=0; % here you can explicitly set the number of left hemisphere clusters  
    default_right_k=0; % number of clusters right hemisphere
    defaultV=5000; 
    % gradientcost, the higher this value is, the higher the enforcement of the gradientprior map will be
    defaultU=1; % datacost, as the costs are relative we can leave this as 1
    % and modify the gradientcost instead. 
    % However, if we want the likelihood of the data to be enforced we can increase this.
    defaultruns=1; % how many  initializations should be started. 
    % For example we are finished with seed 1, should be also compute seed 2, etc.
    defaultexpo=15; % this controls the steepnes of the gradient prior funciton
    % higher values (e.g. 15 or 25) will result in a less stringent gradient map
    % lower values (e.g. 1 or 5) will make more stringend prior map
    defaultiterations=100; % number of optimization iterations 
    default_iter_reduce_gamma=300;
    defaultrandom_init=true; % random initialization, this option might also removed in the future
    defaultgrad_prior='gordon'; % what gradient map we use, default is the map from gordon as in the paper
    defaultfsaverage='fsaverage6'; % which space we use, paper: fsaverage6
    defaultpotts=false; % use no gradient map but potts model, paper setting is false!
    default_lh_avg_file=fullfile(getenv('CBIG_REPDATA_DIR'),'stable_projects','brain_parcellation',...
    'Schaefer2018_LocalGlobal','Time_Matrix_Group_release_full_lh_cov.mat'); 
    %input data, this will only work for CBIG lab
    default_rh_avg_file=fullfile(getenv('CBIG_REPDATA_DIR'),'stable_projects','brain_parcellation',...
    'Schaefer2018_LocalGlobal','Time_Matrix_Group_release_full_rh_cov.mat'); 
    %input data, this will only work for CBIG lab 
    % if you want to compute your own input data, look into CBIG_gwMRF_build_time_matrix and CBIG_gwMRF_build_cov_matrix
    default_initialization_prior='./input/prior_initializations_all.mat'; 
    % this should be given externally with the future prams struct
    default_initialization_motor='./input/clustering/motor_labels.mat';
     % this is not relevant, can be removed
    default_watershed_files='./input/water_files.mat'; 
    % this is only relevant if you want to use a watershed gradient, we did not do this in the paper
    default_output_folder='./output/'; % this is the outputfolder, you may change this everytime
    default_seperate_hemispheres=1; % 0 or 1, we set 1 and compute each hemisphere seperatly
    default_local_concentration=0.0;
    default_estimate_gamma=1;
    default_reduce_gamma=1; 
    default_graphCutIterations=1000; % number of iterations for the graph cut optimizer
    default_first_gamma=5000000; % when there is a split
    default_start_gamma=2*default_first_gamma; % even before a split
    default_kappa_vector=1; % use a different kappa for each different kappa
    default_data_normalization=0; % data normalization is memory intensive
    % if you have large input data sets you should normalize data before hand
    default_left_skip=0; % skip if you want to compute only left hemisphere
    default_right_skip=0; % skip if you want to compute only right hemisphere
    default_reduce_speed=5; % how fast we reduce gamma, with 5 this is gamma'=gamma *1/5
    default_pca=1; % I called this pca, better would be premultiplied matrix
    default_dim=308640; % dimensionality of the input data, relevant for the premultiplied matrix
    % you have to set this, as in the premultiplied matrix there is no way of estimating this from the data
    default_alpha=0.5;
    
    %% These are the key words to assign you own values
    %% assign default values
    addParameter(p,'start_index',defaultseed);
    addParameter(p,'cluster',defaultk);
    addParameter(p,'smoothcost',defaultV);
    addParameter(p,'datacost',defaultU);
    addParameter(p,'runs',defaultruns);
    addParameter(p,'exponential',defaultexpo);
    addParameter(p,'iterations',defaultiterations);
    addParameter(p,'iter_reduce_gamma',default_iter_reduce_gamma);
    addParameter(p,'random_init',defaultrandom_init);
    addParameter(p,'grad_prior',defaultgrad_prior);
    addParameter(p,'fsaverage',defaultfsaverage);
    addParameter(p,'potts',defaultpotts);
    addParameter(p,'lh_avg_file',default_lh_avg_file);
    addParameter(p,'rh_avg_file',default_rh_avg_file);
    addParameter(p,'initialization_prior',default_initialization_prior);
    addParameter(p,'initialization_motor',default_initialization_motor);
    addParameter(p,'output_folder',default_output_folder);
    addParameter(p,'seperate_hemispheres',default_seperate_hemispheres);
    addParameter(p,'local_concentration',default_local_concentration);
    addParameter(p,'left_cluster',default_left_k);
    addParameter(p,'right_cluster',default_right_k);
    addParameter(p,'estimate_gamma',default_estimate_gamma);
    addParameter(p,'graphCutIterations',default_graphCutIterations);
    addParameter(p,'first_gamma',default_first_gamma);
    addParameter(p,'start_gamma',default_start_gamma);
    addParameter(p,'kappa_vector',default_kappa_vector);
    addParameter(p,'data_normalization',default_data_normalization);
    addParameter(p,'reduce_gamma',default_reduce_gamma);
    addParameter(p,'skip_left',default_left_skip);
    addParameter(p,'skip_right',default_right_skip);
    addParameter(p,'reduce_speed',default_reduce_speed);
    addParameter(p,'pca',default_pca);
    addParameter(p,'dim',default_dim);
    addParameter(p,'alpha',default_alpha);
    addParameter(p,'watershed_files',default_watershed_files);

    parse(p,varargin{:})
    prams=p.Results;
    if (prams.local_concentration>0.0)
        default_output_name=['Graph_Cut_faster_','_grad_prior_',prams.grad_prior,'_cluster_',...
        num2str(prams.cluster),'_datacost_',num2str(prams.datacost),'_smoothcost_',...
        num2str(prams.smoothcost),'_iterations_',num2str(prams.iterations),...
        '_local_concentration_',num2str(prams.local_concentration)]
    else
        default_output_name=['Graph_Cut_faster_','_grad_prior_',prams.grad_prior,'_cluster_',...
        num2str(prams.cluster),'_datacost_',num2str(prams.datacost),'_smoothcost_',...
        num2str(prams.smoothcost),'_iterations_',num2str(prams.iterations)]
    end
    addParameter(p,'output_name',default_output_name);
    if(~exist(prams.output_folder)) % create folder if does not exist
        mkdir(prams.output_folder)
    end
    if(~exist([prams.output_folder,'/print/'])) 
        % create folder if does not exist, we will write intermediate results in this
        mkdir([prams.output_folder,'/print/'])
    end
    prams.start_index=check_ischar(prams.start_index);
    prams.runs=check_ischar(prams.runs);

    fileID=fopen([prams.output_folder,'/print/print_',default_output_name,'_seed_',...
    num2str(prams.start_index),'_to_seed_',num2str((prams.start_index+prams.runs-1)),'.txt'],'w');
    %%% open file for output
    [prams.output_folder,'/print/print_',default_output_name,'_seed_',num2str(prams.start_index),...
    '_to_seed_',num2str((prams.start_index+prams.runs-1)),'.txt']
    % Set gradient prior
    % different gradient prioer were tested in the research project.
    % the paper used option 'gordon'
    if(strcmp(prams.grad_prior,'midthickness'))
        fprintf(fileID,'using midthickness gradient prior\n')
        default_lh_grad_file='./lib/input/3_smooth_lh_borders_744_midthickness_subjects_3_postsmoothing_6.mat';
        default_rh_grad_file='./lib/input/3_smooth_rh_borders_744_midthickness_subjects_3_postsmoothing_6.mat';
    elseif(strcmp(prams.grad_prior,'smoothwm'))
        fprintf(fileID,'using smootwhm gradient prior\n')
        default_lh_grad_file='./lib/input/3_smooth_lh_borders_744_smoothwm_subjects_3_postsmoothing_6.mat';
        default_rh_grad_file='./lib/input/3_smooth_rh_borders_744_smoothwm_subjects_3_postsmoothing_6.mat';
    elseif(strcmp(prams.grad_prior,'combined'))
        fprintf(fileID,'using combined gradient prior\n')
        default_lh_grad_file='./lib/input/3_smooth_lh_borders_744_combined_subjects_3_postsmoothing_6.mat';
        default_rh_grad_file='./lib/input/3_smooth_rh_borders_744_combined_subjects_3_postsmoothing_6.mat';
    elseif(strcmp(prams.grad_prior,'gordon'))
        fprintf(fileID,'using gordon gradient prior\n')
        default_lh_grad_file='./lib/input/3_smooth_lh_borders_120_gordon_subjects_3_postsmoothing_6.mat';
        default_rh_grad_file='./lib/input/3_smooth_rh_borders_120_gordon_subjects_3_postsmoothing_6.mat';
    elseif(strcmp(prams.grad_prior,'gordon_min'))
        fprintf(fileID,'using gordon gradient prior\n')
        default_lh_grad_file='./lib/input/3_smooth_lh_borders_120_gordon_subjects_3_postsmoothing_6_min.mat';
        default_rh_grad_file='./lib/input/3_smooth_rh_borders_120_gordon_subjects_3_postsmoothing_6_min.mat';    
    elseif(strcmp(prams.grad_prior,'gordon_fs5'))
        fprintf(fileID,'using gordon gradient prior\n')
        default_lh_grad_file='./lib/input/3_smooth_lh_borders_120_gordon_subjects_3_postsmoothing_6_fs5.mat';
        default_rh_grad_file='./lib/input/3_smooth_rh_borders_120_gordon_subjects_3_postsmoothing_6_fs5.mat';
    elseif(strcmp(prams.grad_prior,'gordon_water'))
        fprintf(fileID,'using gordon water gradient prior\n')
        default_lh_grad_file='./lib/input/3_smooth_lh_borders_120_gordon_subjects_3_postsmoothing_6_water.mat';
        default_rh_grad_file='./lib/input/3_smooth_rh_borders_120_gordon_subjects_3_postsmoothing_6_water.mat';    
    else
        fprintf(fileID,'Warning: unknown gradient prior, will use gordon\n')
        default_lh_grad_file='./lib/input/3_smooth_lh_borders_120_gordon_subjects_3_postsmoothing_6.mat';
        default_rh_grad_file='../lib/input/3_smooth_rh_borders_120_gordon_subjects_3_postsmoothing_6.mat';
    end
    addParameter(p,'lh_grad_file',default_lh_grad_file);
    addParameter(p,'rh_grad_file',default_rh_grad_file);
    parse(p,varargin{:})
    prams=p.Results;
    prams.fileID=fileID;

    %% print some infos
    fprintf(prams.fileID,['outputname:',prams.output_name ,'\n'])
    fprintf(prams.fileID,'using profile %s and %s\n', prams.lh_avg_file, prams.rh_avg_file)
   
    
    if (prams.start_index==defaultseed)
        warning('You use the default seed which makes only sense for debugging')
    end

    %test if inputfiles exist
    if (exist(prams.lh_grad_file,'file')==0)
        error('inputfile %s does not exist \n',prams.lh_grad_file )
    end
    if (exist(prams.rh_grad_file,'file')==0 )
        error('inputfile %s does not exist \n',prams.rh_grad_file )
    end
    if (exist(prams.lh_avg_file,'file')==0)
        error('inputfile %s does not exist \n',prams.lh_avg_file)
    end
    if (exist(prams.rh_avg_file,'file')==0)
        error('inputfile %s does not exist \n',prams.rh_avg_file )
    end


    %% check if input is character if so transform to number
    prams.start_index=check_ischar(prams.start_index);
    prams.runs=check_ischar(prams.runs);
    prams.cluster=check_ischar(prams.cluster);
    prams.smoothcost=check_ischar(prams.smoothcost);
    prams.datacost=check_ischar(prams.datacost);
    prams.exponential=check_ischar(prams.exponential);
    prams.iterations=check_ischar(prams.iterations);
    prams.iter_reduce_gamma=check_ischar(prams.iter_reduce_gamma);
    prams.seperate_hemispheres=check_ischar(prams.seperate_hemispheres);
    prams.local_concentration=check_ischar(prams.local_concentration);
    prams.left_cluster=check_ischar(prams.left_cluster);
    prams.right_cluster=check_ischar(prams.right_cluster);
    prams.start_gamma=check_ischar(prams.start_gamma);
    prams.first_gamma=check_ischar(prams.first_gamma);
    prams.kappa_vector=check_ischar(prams.kappa_vector);
    prams.data_normalization=check_ischar(prams.data_normalization);
    prams.reduce_gamma=check_ischar(prams.reduce_gamma);
    prams.skip_left=check_ischar(prams.skip_left);
    prams.skip_right=check_ischar(prams.skip_right);
    prams.reduce_speed=check_ischar(prams.reduce_speed);
    prams.pca=check_ischar(prams.pca);
    prams.dim=check_ischar(prams.dim);
    prams.alpha=check_ischar(prams.alpha);
    
  
    %%% if left/right cluster are not set, we set it symmetric with value k
    if(prams.left_cluster==0)
       prams.left_cluster=prams.cluster;
       prams.right_cluster=prams.cluster;
    end
    
    %% check if values are in the correct range
    if (prams.local_concentration <0)
        error('local concentration should be non negative')
    end

end

function x=check_ischar(x)
    if(ischar(x))
        x=str2num(x);
    end
end
