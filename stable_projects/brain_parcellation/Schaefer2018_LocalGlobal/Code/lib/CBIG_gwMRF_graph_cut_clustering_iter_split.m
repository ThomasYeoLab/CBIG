function [results]=CBIG_gwMRF_graph_cut_clustering_iter_split(prams)

    % This function loads the data, normalizes it and calls the clustering function
    %
    % Input
    %   - prams = a struct containing various input parameters
    %
    % Ouput
    %   - results = a struct containing the labels and various additional results
    %
    % Example
    %   - [results] = CBIG_gwMRF_graph_cut_clustering_iter_split(prams)
    %
    %Written by Alexander Schaefer and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', prams.fsaverage, 'inflated','cortex');% to identify medial wall/non cortex areas
    l1 = find(lh_avg_mesh.MARS_label == 2);
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', prams.fsaverage, 'inflated', 'cortex');
    r1 = find(rh_avg_mesh.MARS_label == 2);
    prams.dim=prams.dim-1;%to account for lashkari(2010) approximation
    if(strcmp(prams.grad_prior,'gordon_water'))
        CBIG_build_LogOdds(lh_avg_mesh,rh_avg_mesh,prams);
    end

    if((prams.seperate_hemispheres==1)||(prams.local_concentration>0.00))
        %%%% start left hemisphere
        %lh_avg_input = MRIread(prams.lh_avg_file);%% read in the connectivity data
        %lh_vol=reshape(lh_avg_input.vol, lh_num_verts, size(lh_avg_input.vol, 4));
        if(prams.skip_left==0)
            prams.lh_avg_file
            if(prams.pca==0)
                load(prams.lh_avg_file)
                if (size(lh_vol,1)>length(l1))%if rows include the medial wall
                    lh_vol = lh_vol([l1],:);
                end
            else
                lh_vol=load(prams.lh_avg_file);
            end
            clear lh_avg_input;

            rand('twister',5489) % this rand setting is just a double assurance, later the random seed will be defined by prams.seed
            fprintf(prams.fileID,'Computing parameters in split brains \n')
            prams.cluster=prams.left_cluster; %% assign cluster number
            if(prams.pca==1)
                [likeli, results] = CBIG_gwMRF_graph_cut_clustering_split_newkappa_prod(lh_vol,prams,'lh');% start clustering
            else
                [likeli, results] = CBIG_gwMRF_graph_cut_clustering_split_newkappa(lh_vol,prams,'lh');% start clustering
            end
            results.lh_label=results.full_label;
            results.lh_final_likeli=results.final_likeli;
            results.lh_likeli_pos=results.likeli_pos;
            clear results.final_likeli
            clear lh_vol;

        else
            results.lh_label=zeros(size(lh_avg_mesh.MARS_label));
            results.lh_final_likeli=zeros(size(lh_avg_mesh.MARS_label));
            results.lh_likeli_pos=zeros(size(lh_avg_mesh.MARS_label));
            results.D=0;
            results.S=0;
            results.E=0;
            results.UnormalizedE=0;
            results.gamma=zeros(1,prams.left_cluster);
            results.kappa=zeros(1,prams.left_cluster);
            results.mu=0;
        end


        %%%% start right hemisphere
        %rh_avg_input = MRIread(prams.rh_avg_file);
        %rh_vol=reshape(rh_avg_input.vol, rh_num_verts, size(rh_avg_input.vol, 4));
        if(prams.skip_right==0)
            if(prams.pca==0)
                load(prams.rh_avg_file)
                if (size(rh_vol,1)>length(r1))%if rows include the medial wall
                    rh_vol = rh_vol([r1],:);
                end
            else
                rh_vol=load(prams.rh_avg_file);
            end

            clear rh_avg_mesh;% clean some memory
            clear rh_avg_input;
            rand('twister',5489) % this rand setting is just a double assurance, later the random seed will be defined by prams.seed
            prams.cluster=prams.right_cluster;
            if(prams.pca==1)
                [likeli, results_rh] = CBIG_gwMRF_graph_cut_clustering_split_newkappa_prod(rh_vol,prams,'rh');% start clustering
            else
                [likeli, results_rh] = CBIG_gwMRF_graph_cut_clustering_split_newkappa(rh_vol,prams,'rh');% start clustering
            end
            results.rh_label=results_rh.full_label; %% assign cluster number
        else
            results.rh_label=zeros(size(lh_avg_mesh.MARS_label));
            results_rh.final_likeli=zeros(size(lh_avg_mesh.MARS_label));
            results_rh.likeli_pos=zeros(size(lh_avg_mesh.MARS_label));
            results_rh.D=0;
            results_rh.S=0;
            results_rh.E=0;
            results_rh.UnormalizedE=0;
            results_rh.gamma=zeros(1,prams.right_cluster);
            results_rh.kappa=zeros(1,prams.right_cluster);
            results_rh.mu=0;
        end
        clear results.full_label
        results.D=results_rh.D+results.D;
        results.rh_S=results_rh.S;
        results.lh_S=results.S;
        results.S=results_rh.S+results.S;
        results.E=results_rh.E+results.E;
        results.UnormalizedE=results_rh.UnormalizedE+results.UnormalizedE;
        results.gamma=[results.gamma,results_rh.gamma];
        results.kappa=[results.kappa,results_rh.kappa];
        results.rh_final_likeli=results_rh.final_likeli;
        results.rh_likeli_pos=results_rh.likeli_pos;
        clear results.likeli_pos

    end
end

function [sucess]=CBIG_build_LogOdds(lh_avg_mesh6,rh_avg_mesh6,prams)
    load(prams.watershed_files);
    lh_sub_mesh=lh_avg_mesh6;
    rh_sub_mesh=rh_avg_mesh6;
    lh_sub_mesh.MARS_label=(lh_segmentline==-1)+1;
    rh_sub_mesh.MARS_label=(rh_segmentline==-1)+1;
    lh_sub_mesh.MARS_ct.numEntries=2;
    rh_sub_mesh.MARS_ct.numEntries=2;
    [~, ~, lh_prob_mat] = MARS_computeLogOdds(lh_sub_mesh,prams.alpha,Inf);
    [~, ~, rh_prob_mat] = MARS_computeLogOdds(rh_sub_mesh,prams.alpha,Inf);
    [lh_grad_matrix,rh_grad_matrix]=alex_gradient_vertices_to_matrix(lh_prob_mat(2,:),rh_prob_mat(2,:),prams.fsaverage);
    border_matrix=lh_grad_matrix';
    border=mean(lh_grad_matrix',1);
    save(prams.lh_grad_file,'border_matrix','border');

    border_matrix=rh_grad_matrix';
    border=mean(rh_grad_matrix',1);
    save(prams.rh_grad_file,'border_matrix','border');
    sucess=1;
end
