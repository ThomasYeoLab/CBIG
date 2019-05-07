function dxcurr = CBIG_MMLDA_dxchange2dxcurr(dxchange)
% dxcurr = CBIG_MMLDA_dxchange2dxcurr(dxchange)
%
% Convert dxchange variable to dxcurr variable of diagnosis of the subject.
%
% Input:
%   - dxchange  : The changing diagnosis of the patient including 9 cases:
%                 dxchange = 1 means convert from CN to CN
%                 dxchange = 2 means convert from MCI to MCI
%                 dxchange = 3 means convert from AD to AD
%                 dxchange = 4 means convert from CN to MCI
%                 dxchange = 5 means convert from MCI to AD
%                 dxchange = 6 means convert from CN to AD
%                 dxchange = 7 means convert from MCI to CN
%                 dxchange = 8 means convert from AD to MCI
%                 dxchange = 9 means convert from AD to CN
%
% Output:
%   - dxcurr    : The current diagnosis of the patient including 3 cases: 
%                 dxcurr = 1 means Cognitive Normal (CN)
%                 dxcurr = 2 means Mild Cognitive Impairment (MCI)
%                 dxcurr = 3 means Alzheimer's disease (AD)
%
% Example:
%   dxcurr = CBIG_MMLDA_dxchange2dxcurr(4);
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if isnan(dxchange)
    dxcurr = NaN;
else
    switch dxchange
        % Stable: CN to CN
        case 1
            dxcurr = 1;
        % Stable: MCI to MCI
        case 2
            dxcurr = 2;
        % Stable: AD to AD
        case 3
            dxcurr = 3;
        % Conv: CN to MCI
        case 4
            dxcurr = 2;
        % Conv: MCI to AD
        case 5
            dxcurr = 3;
        % Conv: CN to AD
        case 6
            dxcurr = 3;
        % Rev: MCI to CN
        case 7
            dxcurr = 1;
        % Rev: AD to MCI
        case 8
            dxcurr = 2;
        % Rev: AD to CN
        case 9
            dxcurr = 1;
        otherwise
            error('Error: Invalid input dxchange!');
    end
end
