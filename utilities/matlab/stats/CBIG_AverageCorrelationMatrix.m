function CBIG_AverageCorrelationMatrix(SUBJECTS_DIR, input_file, session_file, output_file)

% CBIG_AverageCorrelationMatrix(SUBJECTS_DIR, input_file, session_file, output_file)
% This function is specific to Thomas's project. Please use CBIG_AverageCorrelationMatrixGeneral
% 
% Example:
% CBIG_AverageCorrelationMatrix('/ncf/cnl/20/users/ythomas/MultisessionsSubs/', 'frontal4mm2frontal4mm', 'subjects_session1', 'AvgCorrFrontal4mm2Frontal4mm');
% CBIG_AverageCorrelationMatrix('/ncf/cnl/20/users/ythomas/MultisessionsSubs/', 'frontal4mm2frontal4mm', 'subjects_session2', 'AvgCorrFrontal4mm2Frontal4mm');
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



% Read session 1
fid = fopen([SUBJECTS_DIR '/scripts/' session_file '.txt'], 'r');
count = 0;
while(1)
   tmp = fscanf(fid, '%s\n', 1);
   if(isempty(tmp))
       break
   else
       count = count + 1;
       session{count} = tmp;
   end
end
fclose(fid);

for i = 1:count
    disp([num2str(i) ': ' session{i}]);
    x = load([SUBJECTS_DIR '/' session{i} '/fcMRI_ANALYSIS/' session{i} '_' input_file '.mat']);
    
   if(i == 1)
       corr_mat = x.corr_mat;
   else
       corr_mat = corr_mat + x.corr_mat;
   end
end
corr_mat = corr_mat/count;
save([SUBJECTS_DIR '/' output_file '.' session_file '.mat'], 'corr_mat');
