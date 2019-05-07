function baseline = CBIG_MMLDA_hypotest_dependent_factors(k, N, betas, stats, output, CI_scale)
% baseline = CBIG_MMLDA_hypotest_dependent_factors(k, N, betas, stats, output, CI_scale)
%
% This function will do omnibus and pairwise tests of k factors. 
%
% Input:
%   - k         : number of factors
%   - N         : number of subjects
%   - betas     : beta coefficient from glmfit
%   - stats     : stats from glmfit
%   - output    : txt file
%   - CI_scale  : (optional) confidence interval scale. Default 1.96 (95% CI)
%
% Output:
%   - baseline  : C x 3 matrix. C is number of comparisons between factors.
%                 The 1st column is mean difference between factors.
%                 The 2nd column is stand deviation times the CI_scale
%                 The 3rd column is the p value 
%
% Example:
%   baseline = CBIG_MMLDA_hypotest_dependent_factors(3, 100, betas, stats, './')
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin < 6
    CI_scale = 1.96; % 95% confidence interval
end

num_nuisance = size(betas, 1) - k;
zero_nuisance = zeros(1, num_nuisance);

switch k
    case 2
        % factor2 - factor1
        H = [0 1 zero_nuisance];
        c = 0; % H*b = c

        p = linhyptest(betas, stats.covb, c, H, stats.dfe);

        % output to console and txt file
        fid = fopen(output, 'w');
        fprintf(fid, '%d\n', N)
        if H*betas > 0
            fprintf('factor2 > factor1, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        else
            fprintf('factor2 < factor1, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        end
        fclose(fid);

        idx1 = find(H==1);
        baseline = [H*betas...
            CI_scale*sqrt(stats.covb(idx1, idx1)) p];
    case 3
        % Overall (F2 - F1 = 0; F3 - F1 = 0)
        H = [0 1 0 zero_nuisance; 0 0 1 zero_nuisance];
        c = [0; 0]; % H*b = c
        p = linhyptest(betas, stats.covb, c, H, stats.dfe);

        % output to console and txt file
        fid = fopen(output, 'w')
        fprintf(fid, '%d\n', N)
        fprintf('Overall p = %e\n', p);
        fprintf(fid, '%e\n', p);

        % factor2 - factor1 
        H = [0 1 0 zero_nuisance];
        c = 0; % H*b = c
        p = linhyptest(betas, stats.covb, c, H, stats.dfe);
        if H*betas > 0
            fprintf('factor2 > factor1, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        else
            fprintf('factor2 < factor1, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        end
        idx1 = find(H==1);
        baseline(1, :) = [H*betas...
            CI_scale*sqrt(stats.covb(idx1, idx1))... 
            p];

        % factor3 - factor1
        H = [0 0 1 zero_nuisance];
        c = 0; % H*b = c
        p = linhyptest(betas, stats.covb, c, H, stats.dfe);
        if H*betas > 0
            fprintf('factor3 > factor1, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        else
            fprintf('factor3 < factor1, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        end
        idx1 = find(H==1);
        baseline(2, :) = [H*betas...
            CI_scale*sqrt(stats.covb(idx1, idx1))...
            p];

        % factor3 - factor2
        H = [0 -1 1 zero_nuisance];
        c = 0; % H*b = c
        p = linhyptest(betas, stats.covb, c, H, stats.dfe);
        if H*betas > 0
            fprintf('factor3 > factor2, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        else
            fprintf('factor3 < factor2, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        end
        idx1 = find(H==1);
        idx2 = find(H==-1);
        baseline(3, :) = [H*betas...
            CI_scale*sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)-2*stats.covb(idx1, idx2))...
            p]; 
    case 4
        % Overall (F2 - F1 = 0; F3 - F1 = 0; F4 - F1 = 0)
        H = [0 1 0 0 zero_nuisance; 0 0 1 0 zero_nuisance; 0 0 0 1 zero_nuisance];
        c = [0; 0; 0]; % H*b = c
        p = linhyptest(betas, stats.covb, c, H, stats.dfe);

        % output to console and txt file
        fid = fopen(output, 'w')
        fprintf(fid, '%d\n', N)
        fprintf('Overall p = %e\n', p);
        fprintf(fid, '%e\n', p);


        % factor2 - factor1 
        H = [0 1 0 0 zero_nuisance];
        c = 0; % H*b = c
        p = linhyptest(betas, stats.covb, c, H, stats.dfe);
        if H*betas > 0
            fprintf('factor2 > factor1, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        else
            fprintf('factor2 < factor1, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        end
        idx1 = find(H==1);
        baseline(1, :) = [H*betas...
            CI_scale*sqrt(stats.covb(idx1, idx1))... 
            p];

        % factor3 - factor1
        H = [0 0 1 0 zero_nuisance];
        c = 0; % H*b = c
        p = linhyptest(betas, stats.covb, c, H, stats.dfe);
        if H*betas > 0
            fprintf('factor3 > factor1, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        else
            fprintf('factor3 < factor1, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        end
        idx1 = find(H==1);
        baseline(2, :) = [H*betas...
            CI_scale*sqrt(stats.covb(idx1, idx1))...
            p];

        % factor4 - factor1
        H = [0 0 0 1 zero_nuisance];
        c = 0; % H*b = c
        p = linhyptest(betas, stats.covb, c, H, stats.dfe);
        if H*betas > 0
            fprintf('factor4 > factor1, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        else
            fprintf('factor4 < factor1, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        end
        idx1 = find(H==1);
        baseline(3, :) = [H*betas...
            CI_scale*sqrt(stats.covb(idx1, idx1))...
            p];

        % factor3 - factor2
        H = [0 -1 1 0 zero_nuisance];
        c = 0; % H*b = c
        p = linhyptest(betas, stats.covb, c, H, stats.dfe);
        if H*betas > 0
            fprintf('factor3 > factor2, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        else
            fprintf('factor3 < factor2, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        end
        idx1 = find(H==1);
        idx2 = find(H==-1);
        baseline(4, :) = [H*betas...
            CI_scale*sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)-2*stats.covb(idx1, idx2))...
            p];

        % factor4 - factor2
        H = [0 -1 0 1 zero_nuisance];
        c = 0; % H*b = c
        p = linhyptest(betas, stats.covb, c, H, stats.dfe);
        if H*betas > 0
            fprintf('factor4 > factor2, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        else
            fprintf('factor4 < factor2, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        end
        idx1 = find(H==1);
        idx2 = find(H==-1);
        baseline(5, :) = [H*betas...
            CI_scale*sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)-2*stats.covb(idx1, idx2))...
            p];

        % factor4 - factor3
        H = [0 0 -1 1 zero_nuisance];
        c = 0; % H*b = c
        p = linhyptest(betas, stats.covb, c, H, stats.dfe);
        if H*betas > 0
            fprintf('factor4 > factor3, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        else
            fprintf('factor4 < factor3, p = %e\n', p);
            fprintf(fid, '%e\n', p);
        end
        idx1 = find(H==1);
        idx2 = find(H==-1);
        baseline(6, :) = [H*betas...
            CI_scale*sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)-2*stats.covb(idx1, idx2))...
            p];
            
    otherwise
        warning(['k = ' num2str(k) ' is not specified in current script, you can add it in current script by yourself.'])
end