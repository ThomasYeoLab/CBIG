function [] = CBIG_gwMRF_run_Schaefer2018_create_all_public_parcellations(base_folder)

% [] = CBIG_gwMRF_run_Schaefer2018_create_all_public_parcellations(base_folder)
%
% This script takes all parcellations from Schaefer2018 and maps them from
% fsaverage to different spaces like fslr_32k and MNI space.
% The parcellations in fsaverage space are in the subfolder `input/FreeSurfer5.3/fsaverage/label/`

%Input: 
% - base_folder = String that describes a basefolder that is very specific to
%   the Schaefer2018 project and will contain all output files
% 
%Output:
% - generated as files
%
%Example:
% CBIG_gwMRF_run_Schaefer2018_create_all_public_parcellations([getenv('CBIG_CODE_DIR'),'/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/'])
%
%Written by Alexander Schaefer and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    for i=7:10:17  %% read highres and write low res versions  
        for j=[100,200,400,600,800,1000]
            lh_annot=['input/FreeSurfer5.3/fsaverage/label/lh.Schaefer2018_',num2str(j),'Parcels_',num2str(i) ,'Networks_order.annot'];
            rh_annot=['input/FreeSurfer5.3/fsaverage/label/rh.Schaefer2018_',num2str(j),'Parcels_',num2str(i) ,'Networks_order.annot'];
            copyfile('input/FreeSurfer5.3/',base_folder);
            [a,rh,rh_s]=read_annotation(rh_annot);
            [a,lh,lh_s]=read_annotation(lh_annot);


            for k=2:j/2+1;%%check for usage, we will ignore the first parcel so we start with 2
                if(length(find(lh==lh_s.table(k,5)))<1)
                    warning('parcel in lh not used %i',k)
                end
                if(length(find(rh==rh_s.table(k,5)))<1)
                    warning('parcel in rh not used %i',k)
                end
            end

            mkdir([base_folder,'/FreeSurfer5.3/fsaverage5/label/'])
            mkdir([base_folder,'/FreeSurfer5.3/fsaverage6/label/'])
            mkdir([base_folder,'/HCP/fslr32k/cifti/'])


            write_annotation([base_folder,'/FreeSurfer5.3/fsaverage5/label/lh.Schaefer2018_',num2str(j),'Parcels_',num2str(i) ,'Networks_order.annot'], 0:10241, lh(1:10242), lh_s);
            write_annotation([base_folder,'/FreeSurfer5.3/fsaverage5/label/rh.Schaefer2018_',num2str(j),'Parcels_',num2str(i) ,'Networks_order.annot'], 0:10241, rh(1:10242), rh_s);

            write_annotation([base_folder,'/FreeSurfer5.3/fsaverage6/label/lh.Schaefer2018_',num2str(j),'Parcels_',num2str(i) ,'Networks_order.annot'], 0:40961, lh(1:40962), lh_s);
            write_annotation([base_folder,'/FreeSurfer5.3/fsaverage6/label/rh.Schaefer2018_',num2str(j),'Parcels_',num2str(i) ,'Networks_order.annot'], 0:40961, rh(1:40962), rh_s);

            nparcels=j/2;
            lh_new=zeros(size(lh));
            rh_new=zeros(size(rh));
            for l=1:nparcels
                lh_new(lh==lh_s.table(l+1,5))=l;
                rh_new(rh==rh_s.table(l+1,5))=l;
            end

            system('rm -rf ~/temp123/')
            [lh_fslr32k,rh_fslr32k]=CBIG_project_fsaverage2fsLR(lh_new,rh_new,'fsaverage6','label','~/temp123/','20160827');
            CBIG_gwMRF_write_cifti_from_annot(lh_annot,rh_annot,[base_folder,'/HCP/fslr32k/cifti/Schaefer2018_',num2str(j),'Parcels_',num2str(i) ,'Networks_order'],j/2,lh_fslr32k,rh_fslr32k)



            mkdir([base_folder,'/MNI/'])

            CBIG_project_to_MNI(lh_new,rh_new,[base_folder,'/MNI/Schaefer2018_',num2str(j),'Parcels_',num2str(i) ,'Networks_order']) 
            k=1:j+1;
            cell2csv([base_folder,'/MNI/Schaefer2018_',num2str(j),'Parcels_',num2str(i) ,'Networks_order.txt'],[num2cell(1:2*nparcels)',[lh_s.struct_names(2:end);rh_s.struct_names(2:end)],num2cell([lh_s.table(2:nparcels+1,[1:3]);rh_s.table(2:nparcels+1,[1:3])]),num2cell(zeros(2*nparcels,1))],'\t')
            CBIG_gwMRF_create_FSL_LUT([base_folder,'/MNI/Schaefer2018_',num2str(j),'Parcels_',num2str(i) ,'Networks_order.lut'],[lh_s.table(2:end,[1:3]);rh_s.table(2:end,[1:3])])
        end
    end

    function []=CBIG_project_to_MNI(lh_label,rh_label,filename,mni_mask)
    %%%%%%%%
    % Wrapper funtion that projects from fsaverage to MNI
    % using CBIG_Projectfsaverage2MNI
    %%
    % input label vectors in fsaverage or fsaverge6 space
    % filename is the target ouput filename
    % mni_mask is an optional additional mask
    %%%%
        rh_label(rh_label>0)=rh_label(rh_label>0)+max(lh_label);

        if(size(lh_label, 1) ~= 1)
           lh_label=lh_label';
        end

        if(size(rh_label, 1) ~= 1)
            rh_label=rh_label';
        end

        if (max(size(lh_label))==40962)
            lh_mesh7 = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'sphere', 'cortex');
            rh_mesh7 = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'sphere', 'cortex');
            lh_mesh6 = CBIG_ReadNCAvgMesh('lh', 'fsaverage6', 'sphere', 'cortex');
            rh_mesh6 = CBIG_ReadNCAvgMesh('rh', 'fsaverage6', 'sphere', 'cortex');

            lh_labels7=MARS_NNInterpolate_kdTree(lh_mesh7.vertices,lh_mesh6,lh_label);
            rh_labels7=MARS_NNInterpolate_kdTree(rh_mesh7.vertices,rh_mesh6,rh_label);
        else
            lh_labels7=lh_label;
            rh_labels7=rh_label;
        end
        if (nargin==3)
            output = CBIG_Projectfsaverage2MNI(lh_labels7',rh_labels7');%,mni_mask);
        elseif(nargin==4)
            output = CBIG_Projectfsaverage2MNI(lh_labels7',rh_labels7',mni_mask);
        else
            error('provide correct input arguments')
        end

        MRIwrite(output, [filename,'.nii.gz']);
        MNItemplates_dir = fullfile(getenv('CBIG_CODE_DIR'),'data','templates','volume','FSL5.0.8_MNI_templates');
        system(['mri_vol2vol --targ ', [MNItemplates_dir,'/MNI152_T1_2mm_brain.nii.gz'], ' --regheader --mov ',[filename,'.nii.gz'],' --o ',filename,'_FSLMNI152_2mm.nii.gz']);
        system(['mri_vol2vol --targ ', [MNItemplates_dir,'/MNI152_T1_1mm_brain.nii.gz'], ' --regheader --mov ',[filename,'.nii.gz'],' --o ',filename,'_FSLMNI152_1mm.nii.gz']); 
        system(['flirt -applyisoxfm 3 -interp nearestneighbour -in ',filename,'_FSLMNI152_2mm.nii.gz -ref ', MNItemplates_dir, '/MNI152_T1_3mm_brain.nii.gz -out ',filename, '_3mm.nii.gz'])
        system(['3dcalc -overwrite -a ', filename,'_3mm.nii.gz -expr ''step(a)''  -prefix ',filename,'_3mm_mask.nii.gz']);
        system(['rm -f ', filename,'*.lta']);
        system(['rm -f ', filename,'*.reg']);
        system(['rm -f ', filename,'*mask*']);
        system(['mv ', [filename,'.nii.gz'], ' ', [filename,'_conform.nii.gz']]);

        vol=MRIread([filename,'_FSLMNI152_2mm.nii.gz']);
        [a]=hist(reshape(vol.vol,1,109*91*91),max(rh_label)+1);
        size_of_smallest_parcel=min(a)

    end
end
