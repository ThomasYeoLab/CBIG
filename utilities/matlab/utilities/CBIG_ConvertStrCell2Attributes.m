function att = CBIG_ConvertStrCell2Attributes(str_cell, type, sequence)

% att = CBIG_ConvertStrCell2Attributes(str_cell, type, sequence)
% type can be "unique" or "single"
% if "single", then must specify "sequence_cell", e.g. {'LFT', 'AMB',
% 'RHT'}, and the corresponding att will be a N x 1 vector where LFT is set
% to 0, AMB is set to 1 and RHT is set to 2
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



atts = unique(str_cell);
num_atts = length(atts);

if(strcmp(type, 'single'))
   if(num_atts ~= length(sequence))
      error('Number of unique types not equal length of sequence'); 
   end
    
   if(length(intersect(atts, sequence)) ~= num_atts)
      error('atts and sequence does not match up'); 
   end
end

if(strcmp(type, 'unique'))
   att = zeros(length(str_cell), num_atts);
   for i = 1:length(att)
      found = 0;
      for j = 1:length(atts)
         if(strcmp(str_cell{i}, atts{j}))
            att(i, j) = 1; 
            found = found + 1;
         end
      end
      if(found ~= 1)
         error(['found: ' num2str(found)]); 
      end
   end
   att = att(:, 1:end-1); 
elseif(strcmp(type, 'single'))
   att = zeros(length(str_cell), 1);
   for i = 1:length(att)
      found = 0;
      for j = 1:length(sequence)
         if(strcmp(str_cell{i}, sequence{j}))
            att(i) = j - 1; 
            found = found + 1;
         end
      end
      if(found ~= 1)
         error(['found: ' num2str(found)]); 
      end
   end
else
   error('type not recognized');
end



