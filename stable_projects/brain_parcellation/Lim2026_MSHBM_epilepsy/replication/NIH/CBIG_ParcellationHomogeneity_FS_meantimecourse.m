function  homo_with_weight = CBIG_ParcellationHomogeneity_FS_meantimecourse(lh_labels,rh_labels,mesh,...
    lh_input_filename,rh_input_filename)


% This modified version computes per-parcel homogeneity as the mean of
% correlation between the parcel's average time course and each vertex's
% time course in that parcel (instead of all pairwise vertex-vertex corr).
%
% Written by Mervyn Lim Jun Rui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% ----- shape checks / label stitching (unchanged) -----
if(size(lh_labels,2)~=1)
    lh_labels = lh_labels';
end
if(size(rh_labels,2)~=1)
    rh_labels = rh_labels';
end

if(min(rh_labels(rh_labels~=0)) ~= max(lh_labels) + 1)
    rh_labels(rh_labels~=0) = rh_labels(rh_labels~=0) + max(lh_labels);
end
labels=[lh_labels;rh_labels];

lh_filename=read_sub_list(lh_input_filename);
rh_filename=read_sub_list(rh_input_filename);

lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'inflated', 'cortex');
rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'inflated', 'cortex');
num_subs=size(lh_filename,2);

num_verts = size(lh_avg_mesh.vertices,2)+size(rh_avg_mesh.vertices,2);
max_label = max(labels);

for k=1:num_subs
    fprintf('It is subject %g \n',k);
    count=0;
    lh_curr_filename = textscan(lh_filename{k}, '%s');
    lh_curr_filename = lh_curr_filename{1};
    rh_curr_filename = textscan(rh_filename{k}, '%s');
    rh_curr_filename = rh_curr_filename{1};
    
    num_scans = size(lh_curr_filename,1);
    for i=1:num_scans       
        if (~isnan(lh_curr_filename{i}))
            if(~isempty(lh_curr_filename{i}))
                lh_input=lh_curr_filename{i};
                rh_input=rh_curr_filename{i};
                fprintf('filename: %s \n',lh_input);
                
                [~,lh_vol,~]=read_fmri(lh_input);
                [~,rh_vol,~]=read_fmri(rh_input);
                
                vol=[lh_vol;rh_vol];  % (V x T)
                clear lh_vol rh_vol
                
                all_nan=find(sum(abs(vol),2)==0);
                
                % Will store Fisher-z(vertex-to-parcel-mean corr) per vertex
                homo_full_mat=single(zeros(num_verts,1));
                
                % --- Modified block: compute parcel-mean vs vertex correlations ---
                for c=1:max_label
                    idx = find(labels==c);
                    % remove vertices that are all-zero in this scan
                    idx = setdiff(idx, all_nan);
                    labels_size(k,c)=length(idx);
                    if isempty(idx)
                        continue
                    end
                    
                    % a: (T x Vc) vertex time courses for parcel c
                    a = vol(idx,:)';
                    
                    % Standardize per-vertex time courses (demean, L2 norm)
                    a = bsxfun(@minus, a, mean(a,1));
                    denom = sqrt(sum(a.^2, 1));
                    denom(denom==0) = 1; % guard
                    a = bsxfun(@times, a, 1./denom);
                    
                    % Parcel mean time course (T x 1), standardized same way
                    m = mean(a, 2);
                    m = m - mean(m);
                    m_norm = norm(m);
                    if m_norm == 0
                        corr_vec = zeros(size(a,2),1,'single');
                    else
                        m = m ./ m_norm;
                        % Correlate mean with each vertex: (Vc x 1)
                        corr_vec = single(a' * m); % equals Pearson r due to standardization
                    end
                    
                    % Fisher z for run-averaging stability (use CBIG helper)
                    corr_vec = CBIG_StableAtanh(corr_vec);
                    
                    % Save vertex-aligned z-values
                    homo_full_mat(idx,1) = corr_vec;
                    
                    % Store per-parcel vector (in z-space) for run-averaging below
                    corr_NbyN_mat(c).homo = corr_vec;         % now a VECTOR, not a matrix
                    corr_NbyN_index(c).index = idx;            % keep indices as before
                end
                % -----------------------------------------------------------------
               
                % average across scans (in z-space)
                if(i == 1)
                    homo(:,k) = homo_full_mat;                 % vertex-wise z
                    corr_NbyN_allsub(:,k)=corr_NbyN_mat;       % store vectors in struct
                else
                    homo(:,k) = homo(:,k) + homo_full_mat;
                    for c=1:length(corr_NbyN_mat)
                        % sum z-vectors across runs
                        corr_NbyN_allsub(c,k).homo = corr_NbyN_allsub(c,k).homo + corr_NbyN_mat(c).homo;
                    end
                end
                count=count+1;
            end
        end
    end
    % finalize run-averaging (still z-space here)
    homo(:,k)=homo(:,k)/max(count,1);
    for c=1:length(corr_NbyN_mat)
        corr_NbyN_allsub(c,k).homo = corr_NbyN_allsub(c,k).homo / max(count,1);
    end
end

% ----- Convert to r and compute parcel means per subject (modified) -----
for s=1:num_subs
    for c=1:length(corr_NbyN_mat)
        % Back to r from z for the stored vector
        zvec = corr_NbyN_allsub(c,s).homo;
        rvec = tanh(zvec);
        corr_final_mat(c).homo = rvec; % store r-vector
    end

    %% Compute Homogeneity (modified definition)
    for c=1:length(corr_final_mat)
        rvec = corr_final_mat(c).homo;      % (Vc x 1) correlations of vertices to parcel mean
        if isempty(rvec) || (numel(rvec)==0)
            homo_ci(c,1)=0;
        else
            homo_ci(c,1)=mean(rvec);        % average vertex-to-mean correlation within parcel
        end
    end
    clear corr_final_mat
    homo_with_weight(s,1)=sum(labels_size(s,:)*homo_ci)/sum(labels_size(s,:));
end
end

% ---------------- helpers (unchanged) ----------------
function [fmri, vol, vol_size] = read_fmri(fmri_name)
if (isempty(strfind(fmri_name, '.dtseries.nii')))
    fmri = MRIread(fmri_name);
    vol = single(fmri.vol);
    vol_size = size(vol);
    vol = reshape(vol, prod(vol_size(1:length(vol_size)-1)), vol_size(length(vol_size))); 
    fmri.vol = [];
else
    fmri = ft_read_cifti(fmri_name);
    vol = single(fmri.dtseries);
    vol_size = size(vol);
    fmri.dtseries = [];
end
end

function subj_list = read_sub_list(subject_text_list)
fid = fopen(subject_text_list, 'r');
i = 0;
while(1)
    tmp = fgetl(fid);
    if(tmp == -1)
        break
    else
        i = i + 1;
        subj_list{i} = tmp;
    end
end
fclose(fid);
end