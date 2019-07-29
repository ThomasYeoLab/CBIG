function homo_with_weight = CBIG_ComputeParcellationHomogeneity_fslr(lh_labels,rh_labels,input_filename)

% homo_with_weight = CBIG_ComputeParcellationHomogeneity_fslr(lh_labels,rh_labels,input_filename)
%
% This function is used to compute the homogeneity metric with a given
% parcellation lh_labels, rh_labels in fslr space. To compute the
% homogeneity metric, the function will first compute the averaged pairwise 
% correlation c_i within each parcel/network i, assuming the size of the
% network/parcel is s_i and the total number of vertcies is N. The
% homogeneity metric is defined as: sum(c_i x s_i)/N
%
% INPUT:
%       - lh_labels, rh_labels: (N: #vertices)
%         A Nx1 or 1xN vector. Parcellation labels. The medial wall area
%         should be set as label 0.
%
%       - input_filename:
%         A text file. The list of fMRI data with multiple rows, where each
%         row corresponds to a subject, each row contains fMRI data file
%         names of multiple runs. Each run is separated by a white space.
%         For example, the text file can be:
%         subject1_rs_fMRI_run1.nii.gz subject1_rs_fMRI_run2.nii.gz
%         subject2_rs_fMRI_run1.nii.gz
%         subject3_rs_fMRI_run1.nii.gz subject3_rs_fMRI_run2.nii.gz
%
% OUTPUT:
%       - homo_with_weight:
%         A #subject x 1 vector, the homogeneity metric for each subject.
%
% Example:
%
% input_filename = '/data/users/rkong/storage/ruby/data/Individual_areal_HCP/data/split_data/HCP_dataC/lists/fMRI_list_row/102816_rfMRI_postproc_list.txt';
% load('/data/users/rkong/storage/ruby/data/Individual_areal_HCP/split_data/test/HCP/40/none/400/V2_4session/V2a/parcellation/Individual_Intra_MRF_sub1_w20_MRF50.mat');
% homo_with_weight = CBIG_ComputeParcellationHomogeneity_fslr(lh_labels,rh_labels,input_filename)
%
%Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% This function does not need vector check because the function itself
% contains checking statement.

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

filename=read_sub_list(input_filename);

lh_avg_mesh=CBIG_read_fslr_surface('lh','fs_LR_32k','inflated');
rh_avg_mesh=CBIG_read_fslr_surface('rh','fs_LR_32k','inflated');
num_subs=size(filename,2);

for k=1:num_subs
    fprintf('It is subject %g \n',k);
    count=0;
    curr_filename = textscan(filename{k}, '%s');
    curr_filename = curr_filename{1}; % a cell of size (#runs x 1) for subject k in the first list
    num_scans = size(curr_filename,1);
    for i=1:num_scans       
        if (~isnan(curr_filename{i}))
            if(~isempty(curr_filename{i}))
                input=curr_filename{i};
                fprintf('filename: %s \n',input);
                [~, vol, ~] = read_fmri(input);
                vol(64985:end,:)=[];
                
                all_nan=find(isnan(mean(vol,2))==1);
                %left hemi
                homo_full_mat=single(zeros(size(lh_avg_mesh.vertices,2)+size(rh_avg_mesh.vertices,2),1));
                for c=1:max(labels)      
                    a=vol(labels==c,:)';
                    index_cluster=find(labels==c);
                    index_nan=transpose(find(isnan(mean(a))==1));
                    a(:,index_nan)=[];
                    index_cluster=setdiff(index_cluster,all_nan);
                    labels_size(k,c)=length(index_cluster);
                    
                    a=bsxfun(@minus,a,mean(a,1));
                    a=bsxfun(@times,a,1./sqrt(sum(a.^2, 1)));
                    corr_mat=a'*a;
                    corr_mat=CBIG_StableAtanh(corr_mat);

                    %NxN corr matrix, N: # of vertcies within cluster c
                    corr_NbyN_mat(c).homo=corr_mat;
                    %index of NxN matrix
                    corr_NbyN_index(c).index=index_cluster;                 
                    %average of NxN matrix, mean_corr_NbyN_mat:Nx1;
                    mean_corr_NbyN_mat=mean(corr_mat,2);
                    homo_full_mat(index_cluster,1)=mean_corr_NbyN_mat;                    
                end
  
                 %average across scans
                 if(i == 1)
                     homo(:,k) = homo_full_mat;
                     corr_NbyN_allsub(:,k)=corr_NbyN_mat;
                 else
                     homo(:,k) = homo(:,k) + homo_full_mat;
                     for c=1:length(corr_NbyN_mat)
                         corr_NbyN_allsub(c,k).homo=corr_NbyN_allsub(c,k).homo+corr_NbyN_mat(c).homo;
                     end                   
                end
                count=count+1;
            end
        end
    end
    homo(:,k)=homo(:,k)/count;
    for c=1:length(corr_NbyN_mat)
        corr_NbyN_allsub(c,k).homo=corr_NbyN_allsub(c,k).homo/count;
    end
end

for s=1:num_subs
    for c=1:length(corr_NbyN_mat)
        corr_final_mat(c).homo=corr_NbyN_allsub(c,s).homo;
        corr_final_mat(c).homo=tanh(corr_final_mat(c).homo);
    end
    
    %% Compute Homogeneity
    for c=1:length(corr_final_mat)
        temp=corr_final_mat(c).homo;
        homo_ci(c,1)=(sum(sum(temp))-size(temp,1))/(size(temp,1)*(size(temp,1)-1));
        if(size(temp,1)==1||size(temp,1)==0)
            homo_ci(c,1)=0;
        end
    end
    clear corr_final_mat
    homo_with_weight(s,1)=sum(labels_size(s,:)*homo_ci)/sum(labels_size(s,:));
end
end

function [fmri, vol, vol_size] = read_fmri(fmri_name)

% [fmri, vol] = read_fmri(fmri_name)
% Given the name of functional MRI file (fmri_name), this function read in
% the fmri structure and the content of signals (vol).
% 
% Input:
%     - fmri_name:
%       The full path of input file name.
%
% Output:
%     - fmri:
%       The structure read in by MRIread() or ft_read_cifti(). To save
%       the memory, fmri.vol (for NIFTI) or fmri.dtseries (for CIFTI) is
%       set to be empty after it is transfered to "vol".
%
%     - vol:
%       A num_voxels x num_timepoints matrix which is the content of
%       fmri.vol (for NIFTI) or fmri.dtseries (for CIFTI) after reshape.
%
%     - vol_size:
%       The size of fmri.vol (NIFTI) or fmri.dtseries (CIFTI).

if (isempty(strfind(fmri_name, '.dtseries.nii')))
    % if input file is NIFTI file
    fmri = MRIread(fmri_name);
    vol = single(fmri.vol);
    vol_size = size(vol);
    vol = reshape(vol, prod(vol_size(1:3)), prod(vol_size)/prod(vol_size(1:3)));
    fmri.vol = [];
else
    % if input file is CIFTI file
    fmri = ft_read_cifti(fmri_name);
    vol = single(fmri.dtseries);
    vol_size = size(vol);
    fmri.dtseries = [];
end


end

function subj_list = read_sub_list(subject_text_list)
% this function will output a 1xN cell where N is the number of
% subjects in the text_list, each subject will be represented by one
% line in the text file
% NOTE: multiple runs of the same subject will still stay on the same
% line
% Each cell will contain the location of the subject, for e.g.
% '<full_path>/subject1_run1_bold.nii.gz <full_path>/subject1_run2_bold.nii.gz'
fid = fopen(subject_text_list, 'r');
i = 0;
while(1);
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