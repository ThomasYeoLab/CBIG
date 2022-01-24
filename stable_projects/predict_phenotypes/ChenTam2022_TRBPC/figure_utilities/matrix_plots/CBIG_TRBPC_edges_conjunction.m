function CBIG_TRBPC_edges_conjunction(path_data, lim_maps)
% CBIG_TRBPC_edges_conjunction(path_data, lim_maps)
%
% This function will do a conjunction analysis and generate conjunction
% matrices, masked matrices after applying a logical conjunction mask
%
% required input: 
% - path_data : a filepath to a directory that contains a .mat file called
%           `matrices_cat_masked_sig.mat` which has a structure named
%           `mean_matrix_cat_sig` ; the outputs will also be written to this folder
% - lim_maps : a vector setting the lower and upper limits of the output
%           matrices. If set to [], the default values are [-2 2]
%
% Written by Angela Tam & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%% set variables

% path to relevance maps that have been masked out for significant edges
% after permutation testing (contains a structure MEAN_MATRIX_CAT_SIG)
path_mat_sig = [path_data filesep 'matrices_cat_masked_sig.mat'];

% set path to save outputs
path_out = [path_data filesep 'conjunction'];
if ~exist(path_out), mkdir(path_out); end

path_out_chord = [path_data filesep 'conjunction/chord'];
if ~exist(path_out_chord), mkdir(path_out_chord); end

path_out_surface = [path_data filesep 'conjunction/surfaces'];
if ~exist(path_out_surface), mkdir(path_out_surface); end

% set limits for plots of the 419x419 matrices
lims_conj = [3 4]; % limit for conjunction plots
if isempty(lim_maps)
    lim_maps = [-2 2]; % limits for the relevance maps
end

% names of behavioral cluster inside path_mat
clus = {'cog','ment_health','pers'};
% names of fmri modalities inside path_mat
mods = {'rs','mid','nback','sst'};

% load data
load(path_mat_sig)

% generate network names and pairs
[networks,blocks,net_ind] = CBIG_TRBPC_generate_block_pairs_18networks;

%% create logical masks for each behavioral cluster & fMRI condition map

% structure to store results
mean_matrix_logical = [];

for cc = 1:length(clus)
    clu = clus{cc};
    for mm = 1:length(mods)
        mod = mods{mm};
        mat = mean_matrix_cat_sig_std.(clu).(mod);
        %% case 1: absolute
        % logical matrix of the absolute values
        mean_matrix_logical.abs.(clu).(mod) = logical(abs(mat));
        %% case 2 & 3: positive and negative values
        % create empty matrices
        mean_matrix_logical.pos.(clu).(mod) = zeros(size(mat));
        mean_matrix_logical.neg.(clu).(mod) = zeros(size(mat));
        % avg within network block to determine if a block is overall
        % positive or negative
        for bb = 1:length(blocks)
            blk = blocks{bb};
            tmp = strsplit(blk,'_');
            net_a = tmp{1};
            net_b = tmp{2};
            ind_net_a = find(contains(networks,net_a));
            ind_net_b = find(contains(networks,net_b));
            % take that block from the matrix
            corr_blk = mat(net_ind.(net_b),net_ind.(net_a));
            
            % if that block is on the diagonal
            if isequal(net_a,net_b)
                % make a mask of the lower triangle of that block
                mask_blk_tril = tril(ones(size(corr_blk,1),size(corr_blk,2)),-1);
                % average the values of the lower triangle of that block
                corr_vec_tril = corr_blk(mask_blk_tril==1);
                avg_corr = mean(corr_vec_tril);
            else
                % average acorss the entire block
                avg_corr = mean(corr_blk(:));
            end
            
            if avg_corr > 0
                % if avg of network block is positive, fill that block with ones
                mean_matrix_logical.pos.(clu).(mod)(net_ind.(net_b),net_ind.(net_a)) = 1;
                mean_matrix_logical.pos.(clu).(mod)(net_ind.(net_a),net_ind.(net_b)) = 1;
            elseif avg_corr < 0
                % if avg of network block is negative, fill that block with ones
                mean_matrix_logical.neg.(clu).(mod)(net_ind.(net_b),net_ind.(net_a)) = 1;
                mean_matrix_logical.neg.(clu).(mod)(net_ind.(net_a),net_ind.(net_b)) = 1;
            else
                % if neither positive or negative, do nothing
                continue
            end
        end
        % make the positive and negative masks logical
        mean_matrix_logical.pos.(clu).(mod) = logical(mean_matrix_logical.pos.(clu).(mod));
        mean_matrix_logical.neg.(clu).(mod) = logical(mean_matrix_logical.neg.(clu).(mod));
    end  
end

%% do the conjunction analysis

% structures to store results
mat_logical_conjunction_clus = [];
mask_logical_clus_all_fmri = [];

vals = {'pos','neg','abs'};
for vv = 1:length(vals)
    val = vals{vv};
    % sum the logical matrices within the behavioural clusters
    for cc = 1:length(clus)
        clu = clus{cc};
        tmp_mat = zeros(419,419,4);
        % concatenate the matrices for each fmri modality
        for mm = 1:length(mods)
            mod = mods{mm};
            tmp_mat(:,:,mm) = mean_matrix_logical.(val).(clu).(mod);
        end
        mat_logical_conjunction_clus.(val).(clu) = sum(tmp_mat,3);
        % plot the conjunction matrix
        CBIG_TRBPC_plot_fc_matrix(mat_logical_conjunction_clus.(val).(clu),lims_conj,'hot')
        fname = strcat(path_out, filesep, 'conjunction_', val, '_', clu,...
            '_min', num2str(lims_conj(1)), 'max', num2str(lims_conj(2)));
        saveas(gcf,fname,'svg')
        % save mask for connections that are in all four fMRI states
        mask_logical_clus_all_fmri.(val).(clu) = mat_logical_conjunction_clus.(val).(clu) == 4;
    end
end

mat_out = [path_out filesep 'mask_logical_clus_all_fmri.mat'];
save(mat_out, 'mask_logical_clus_all_fmri');

%% average the relevance matrices for each behavioral cluster and apply conjunction masks
% (i.e. use the binary mask from the conjunction to overlay on top of the avg matrices)

vec_relevance_conjunction = [];
masked_avg_relevance_conjunction = [];

% for each case (absolute, positive, negative)
for vv = 1:length(vals)
    val = vals{vv};
    % for each behavioral cluster
    for cc = 1:length(clus)
        clu = clus{cc};
        % get the mask
        mask_clus = mask_logical_clus_all_fmri.(val).(clu);
        %% average the behavioural clusters maps across the 4 fmri conditions
        stacked_mat_tmp = zeros(419,419,4);
        nn = 1;
        % for each fmri modality
        for mm = 1:length(mods)
            mod = mods{mm};
            if strcmp(val, 'abs')
                stacked_mat_tmp(:,:,nn) = abs(mean_matrix_cat_sig_std.(clu).(mod));
                stacked_mat_tmp_raw(:,:,nn) = abs(mean_matrix_cat_sig_raw.(clu).(mod));
            else
                stacked_mat_tmp(:,:,nn) = mean_matrix_cat_sig_std.(clu).(mod);
                stacked_mat_tmp_raw(:,:,nn) = mean_matrix_cat_sig_raw.(clu).(mod);
            end
            nn = nn + 1;
        end
        mat_avg_clus = mean(stacked_mat_tmp,3);
        mat_avg_clus_raw = mean(stacked_mat_tmp_raw,3);
        
        %% apply the conjunction mask
        mat_masked_avg_clus = mat_avg_clus.*mask_clus;
        mat_masked_avg_clus_raw = mat_avg_clus_raw.*mask_clus;
        % plot the matrix that is divided by std for visualization
        CBIG_TRBPC_plot_fc_matrix(mat_masked_avg_clus,lim_maps,'hot_cold')
        fname = strcat(path_out, filesep, 'relevance_masked_conjunction_', val, '_', clu);
        saveas(gcf,fname,'svg')
        % save the masked matrices
        masked_avg_relevance_conjunction.avg.(val).(clu) = mat_masked_avg_clus;
        
        %% average the values within network blocks to prep for chord diagrams
        avg_rel_blk = CBIG_TRBPC_avg_network_pairs_relevance(mat_masked_avg_clus);
        % write to csv
        csv_out = [path_out filesep 'chord/avg_' val '_' clu '.csv'];
        csvwrite(csv_out, avg_rel_blk)
        
        %% grab only the 400 cortical parcels from the average masked matrix and create files to project to surface
        % rearrange the matrix back to original ordering (so it will be
        % compatible with function to project to surface)
        [old_ind,~] = CBIG_TRBPC_newmatrix_2_oldmatrix_order(17);
        mat_masked_avg_clus_o = mat_masked_avg_clus_raw(old_ind,old_ind);
        mat_conj_tmp = mat_masked_avg_clus_o(:,1:400);
        % make absolute
        mat_conj_tmp = abs(mat_conj_tmp);
        % sum the values for each parcel & vectorize
        vec = sum(mat_conj_tmp,1);
        % save this vector into structure
        vec_relevance_conjunction.(val).(clu) = vec;
        
        %% grab the subcortical
        mat_subcort = mat_masked_avg_clus_o(:,401:419);
        % make absolute
        mat_subcort = abs(mat_subcort);
        % sum the values for each parcel & vectorize
        vec_subcort = sum(mat_subcort,1);
        % save this vector
        vec_subcortical.(val).(clu) = vec_subcort;
    end
end

% save the masked matrices
mat_out = [path_out filesep 'mat_relevance_conjunction_masked.mat'];
save(mat_out, 'masked_avg_relevance_conjunction')

% save the vector for projection to surface for cortical ROIs
vec_out = [path_out filesep 'surfaces/vec_relevance_sum_conjunction.mat'];
save(vec_out,'vec_relevance_conjunction')

% save the vector for the subcortical ROIs
vec_s_out = [path_out filesep 'surfaces/vec_relevance_sum_conjunction_subcortical.mat'];
save(vec_s_out,'vec_subcortical')




