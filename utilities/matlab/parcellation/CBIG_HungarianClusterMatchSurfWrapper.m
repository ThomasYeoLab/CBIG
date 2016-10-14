function CBIG_HungarianClusterMatchSurfWrapper(ref_file, input_file, output_file)

% CBIG_HungarianClusterMatchSurfWrapper(ref_file, input_file, output_file)
%
% This is the wrapper function of surface labels Hungarian matching.
%
% Given a reference file name, this function match labels of input file
% with the reference file.
%
% Reference file and input file are .mat file on surface, and they should  
% have "lh_labels" and "rh_labels" fields.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



ref_labels = load(ref_file);
input_labels = load(input_file);

[labels, assign, cost, dice_overlap] = CBIG_HungarianClusterMatch([ref_labels.lh_labels; ref_labels.rh_labels], [input_labels.lh_labels; input_labels.rh_labels]);
disp(num2str(cost));
disp(num2str(dice_overlap));

lh_labels = labels(1:length(input_labels.lh_labels));
rh_labels = labels(length(input_labels.lh_labels)+1:end);

if(isfield(input_labels, 'lh_s'))
    lh_s = input_labels.lh_s;
    rh_s = input_labels.rh_s;
else
    lh_s = [];
    rh_s = [];
end

if(isfield(ref_labels, 'colors'))
    colors = ref_labels.colors;
else
    colors = [];
end

if(isfield(input_labels, 'mtc'))
  mtc = zeros(size(input_labels.mtc));
  p = zeros(size(input_labels.p));
  for i = 1:size(mtc, 1)
      mtc(assign(i), :) = input_labels.mtc(i, :);
      p(assign(i)) = input_labels.p(i);
  end
  lambda = input_labels.lambda;
  save(output_file, 'lh_labels', 'rh_labels', 'cost', 'dice_overlap', 'lambda', 'mtc', 'p', 'lh_s', 'rh_s', 'colors');
else 
  save(output_file, 'lh_labels', 'rh_labels', 'cost', 'dice_overlap', 'lh_s', 'rh_s', 'colors');
end
% exit
