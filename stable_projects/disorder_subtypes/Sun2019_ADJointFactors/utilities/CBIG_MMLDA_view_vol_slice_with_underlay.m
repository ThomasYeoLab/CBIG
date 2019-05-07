function CBIG_MMLDA_view_vol_slice_with_underlay(in_vol, underlay_vol, mni_space, color_map_name, plane, min_thresh, max_thresh, out_dir, out_name)
% CBIG_MMLDA_view_vol_slice_with_underlay(in_vol, underlay_vol, ...
%   mni_space, color_map_name, plane, min_thresh, max_thresh, out_dir, out_name)
%
% Overlay input volume on underlay volumn (e.g., MNI template) with colormap.
% Then, take screeshots of different slices of the volume and save it as images. 
%
% Input:
%   - in_vol        : full path of input volume
%   - underlay_vol  : full path of underlay volumn, e.g., MNI template
%   - mni_space     : "MNI2mm" or "MNI1.5mm"
%   - color_map_name: "MingHot" or "ColdHot". See colorbar in "colormaps" folder
%   - plane         : "coronal" or "sagittal" or "axial" plane
%   - min_thresh    : if "MingHot" minimum of the colorscale
%                     if "ColdHot" absolute minimum of the colorscale 
%   - max_thresh    : if "MingHot" maximum of the colorscale 
%                     if "ColdHot" absolute maximum of the colorscale
%   - out_dir       : output directory of slices 
%   - out_name      : output name of slices
%
% Examples:
% CBIG_MMLDA_view_vol_slice_with_underlay('topic1.nii.gz', 'MNI_Template.nii.gz', ...
%   'MNI2mm', 'MingHot', 'coronal', 7.5e-6, 1.5e-5, '~/storage', 'factor1')
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

mri = MRIread(in_vol);
vols = mri.vol;

switch color_map_name
    case 'MingHot'
        no_bins = 255*3;
        color_map = ['lut:lut=' CBIG_CODE_DIR ...
        '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities/' ...
        'colormaps/MingHotColormap/MingHotColormap.txt'];
        % discretize the input volume into categorical numbers corresponding to different
        % bins in the color map
        edges = min_thresh:((max_thresh-min_thresh)/no_bins):max_thresh;
        ind = vols>min_thresh&vols<max_thresh;
        vols(vols<min_thresh) = 0; % do not display <MIN
        vols(vols>max_thresh) = no_bins; % display >MAX as MAX
        vols(ind) = discretize(vols(ind), edges); % discretize others
    case 'Jet'
        no_bins = 64;
        color_map = ['lut:lut=' CBIG_CODE_DIR ...
        '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities/' ...
        'colormaps/JetColormap/JetColormap.txt'];
        % Discretize to use JetColormap
        ind = abs(vols)<min_thresh;
        vols(vols>max_thresh) = max_thresh;
        vols(vols<-max_thresh) = -max_thresh;
        edges = -max_thresh:(2*max_thresh/no_bins):max_thresh;
        vols = discretize(vols, edges);
        % Up to this point, [-minAbsVal, minAbsVal] is mapped to the middle of the colormap (black color)
        % We don't want to display them, so set them to 0 (min lable is 1)
        vols(ind) = 0;
    case 'ColdHot'
        no_bins = 255*4;
        color_map = ['lut:lut=' CBIG_CODE_DIR ...
        '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities/' ...
        'colormaps/ColdHotColormap/ColdHotColormap.txt'];
        % Discretize to use mingColdHotColormap
        ind = abs(vols)<min_thresh;
        vols(vols>max_thresh) = max_thresh;
        vols(vols<-max_thresh) = -max_thresh;
        edges = -max_thresh:(2*max_thresh/no_bins):max_thresh;
        vols = discretize(vols, edges);
        % Up to this point, [-minAbsVal, minAbsVal] is mapped to the middle of the colormap (black color)
        % We don't want to display them, so set them to 0 (min lable is 1)
        vols(ind) = 0;
    otherwise
        error('No such option.')
end

% write out the discretized volume 
mri.vol = vols;
mkdir(out_dir)
MRIwrite(mri, [out_dir '/' out_name '_discrete.nii.gz']);
overlay_vol = [out_dir '/' out_name '_discrete.nii.gz'];

% take screenshort of volume slice
cmd = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities/CBIG_MMLDA_view_vol_slice_with_underlay.sh ' ...
    overlay_vol ' ' underlay_vol ' ' mni_space ' ' color_map ' ' plane ' ' out_dir ' ' out_name];
disp(cmd)
system(cmd)

