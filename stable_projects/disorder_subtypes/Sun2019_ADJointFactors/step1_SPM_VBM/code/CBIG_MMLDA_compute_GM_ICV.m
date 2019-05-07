% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


addpath(spm_dir)

% Read rid.txt file
fid=fopen(id_list);
rid=textscan(fid,'%s');
fclose(fid);
rid=rid{1,1};

% Calculate the GM and ICV of all subjects
for i = 1:length(rid)
    gm_file  = [output_dir '/mri/p1' rid{i} '.nii'];
    wm_file  = [output_dir '/mri/p2' rid{i} '.nii'];
    csf_file = [output_dir '/mri/p3' rid{i} '.nii'];

    gm_vol(i) = CBIG_MMLDA_single_image_vol(gm_file);
    wm_vol(i) = CBIG_MMLDA_single_image_vol(wm_file);
    csf_vol(i)= CBIG_MMLDA_single_image_vol(csf_file);
end
icv_vol=sum([gm_vol; wm_vol; csf_vol]);
id_gmVol_icv=[rid num2cell(gm_vol') num2cell(icv_vol')];
mkdir([output_dir '/gmVolAndICV'])
save([output_dir '/gmVolAndICV/id_gmVol_icv.mat'], 'id_gmVol_icv');
cell2csv([output_dir '/gmVolAndICV/id_gmVol_icv.csv'], id_gmVol_icv)

rmpath(spm_dir)