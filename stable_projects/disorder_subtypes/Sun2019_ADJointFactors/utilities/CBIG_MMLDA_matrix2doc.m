function  CBIG_MMLDA_matrix2doc(output_filename, input_mtx)
% CBIG_MMLDA_matrix2doc(output_filename, input_mtx)
%
% This function convert input matrix to output document which is an input of
% CBIG_MMLDA_est.sh or CBIG_MMLDA_inf.sh. For the docmument format, please refer to
% $CBIG_CODE_DIR/external_packages/mmlda-c-dist/README.md
%
% Input:
%       out_name: the full path of output doc (e.g., brain.dat)
%       input_mtx: input N x V matrix, N is the # of subjects and V is # of voxels  
% 
% Example:
%       CBIG_MMLDA_matrix2doc('xxx.dat', input_mtx)
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

fileID = fopen(output_filename, 'W'); % open file and write

noSubjects = size(input_mtx, 1);

for idx1 = 1:noSubjects    
    oneSub = input_mtx(idx1, :);
    noTerms = sum(oneSub~=0);
    % print # of words in this documnet
    fprintf(fileID, '%i ', noTerms);
    for idx2 = 1:numel(oneSub)
        % print words whose word counts are not zero
        if oneSub(idx2) ~= 0
            fprintf(fileID, '%i:%i ', idx2-1, oneSub(idx2));
        end
    end
    fprintf(fileID, '\n');
    disp(idx1);        
end
fclose(fileID);