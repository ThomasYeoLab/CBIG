function [map, lab_final] = CBIG_MakeMaximalColormap(n)

% [map, lab_final] = CBIG_MakeMaximalColormap(n)
%
% This function is used to create colormap with n different colors. Colors
% are as different as possible.
%
% Input:
%      -n: (double or int)
%       number of different colors.
%
% Output:
%       -map:
%        nx3 matrix. n is the number of different colors. each row
%        corresponds to the rgb value of the color.
%
%       -lab_final:
%        nx3 matrix. n is the number of different colors. each row
%        corresponds to the lab value of the color. lab_final is used to
%        generate map using map = applycform(lab_final, 'lab2srgb')
%
% Example:
% [map, lab_final] = CBIG_MakeMaximalColormap(8)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(n <= 2)
   error('n should be greater than 2'); 
end

num_split = ceil(n^(1/3));
l = linspace(0, 100, num_split);
a = linspace(-100, 100, num_split);
b = linspace(-100, 100, num_split);

lab = [];
for k = 1:length(b)
    for j = 1:length(a)
        for i = 1:length(l)
            lab = [lab; [l(i) a(j) b(k)]];
        end
    end
end

% Now sort it so that each subsequent color is as far from the rest of the
% color as possible.
lab_final = [0 0 0];
selected_bool = zeros(size(lab, 1), 1);
tmp = find(sum(abs(lab), 2) == 0);
if(~isempty(tmp))
   selected_bool(tmp) = 1; 
end

for i = 2:n
    for j = 1:size(lab, 1)
        dist(j) = sum(sum(abs(repmat(lab(j, :), size(lab_final, 1), 1) - lab_final)));  
    end
    
    [Y, index] = sort(dist, 'descend');
    for j = 1:size(lab, 1)
        if(selected_bool(index(j)) == 0)
            selected_bool(index(j)) = 1;
            lab_final = [lab_final; lab(index(j), :)];
            break;
        end
    end
end



cform = makecform('lab2srgb');
map = applycform(lab_final, cform);