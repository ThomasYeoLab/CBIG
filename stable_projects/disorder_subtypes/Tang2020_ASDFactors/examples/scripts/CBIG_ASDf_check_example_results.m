function CBIG_ASDf_check_example_results(output_dir, log_file)
% CBIG_ASDf_check_example_results(output_dir, log_file)
%
% This function compares the output after you run the example code with the
% reference output to ensure that you have properly downloaded codes in
% Tang2020_ASDFactors and the example code works fine. The comparison results
% will be saved in the log file..
%
% '[FAILED]' message indicates there is something wrong with the code.
% Your output is the same as our reference result if and only if all messages show '[PASSED]'.
%
% Input:
%     - output_dir:
%           Absolute path to the directory where your output of example
%           code is located
%     - log_file:
%           Absolute path to the log file where comparison results will be
%           printed
%
% Example:
%       CBIG_ASDf_check_example_results('~/Tang2019_ASDFactors/example_output',
%          '~/Tang2020_ASDFactors/example_log.txt');
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
               'disorder_subtypes', 'Tang2020_ASDFactors', ...
               'examples', 'results');

%% Open log file
log_fileID = fopen(log_file, 'a');

%% Check z-normalized, discretized FC data
%%% For estimate
% fprintf('\n--------[CHECK 1] Checking FC2doc results (factor estimate):\n');
fprintf(log_fileID,'\n--------[CHECK 1] Checking FC2doc results (factor estimate):\n');

% Check if dx1.dat and dx2.dat exist
if (~exist(fullfile(output_dir,'FC2doc','step1_output_dx1.dat'), 'file'))
%     fprintf('[FAILED] Word counts document for ASD participants missing!\n');
    fprintf(log_fileID, '[FAILED] Word counts document for ASD participants missing!\n');
    error('[FAILED] Word counts document for ASD participants missing!');
end

if (~exist(fullfile(output_dir,'FC2doc','step1_output_dx2.dat'), 'file'))
%     fprintf('[FAILED] Word counts document for control participants missing!\n');
    fprintf(log_fileID, '[FAILED] Word counts document for control participants missing!\n');
    error('[FAILED] Word counts document for control participants missing!');
end

% Check discretized z scores, these are what were written into .dat word
% counts documents.
if (~exist(fullfile(output_dir,'FC2doc','step1_output_zScores.mat'), 'file'))
%     fprintf('[FAILED] Z scores missing!');
    fprintf(log_fileID, '[FAILED] Z scores missing!');
    error('[FAILED] Z scores missing!');
else
    zScores = load(fullfile(output_dir,'FC2doc','step1_output_zScores.mat'));
    ref_zScores = load(fullfile(ref_dir,'FC2doc','step1_output_zScores.mat'));
    
    if (~isequal(ref_zScores.discretized_z, zScores.discretized_z))
        diff_discretZ = max(max(abs(ref_zScores.discretized_z - zScores.discretized_z)));
%         fprintf('[FAILED] Dicretized Z-scores are different from reference, max abs diff = %f.\n]', diff_discretZ);
        fprintf(log_fileID, '[FAILED] Dicretized Z-scores differ from ref, max abs diff=%f.\n]', diff_discretZ);
        error('[FAILED] Dicretized Z-scores differ from ref, max abs diff=%f.]', diff_discretZ);
    else
%         fprintf('[PASSED] Dicretized Z-scores are the same as reference.\n');
        fprintf(log_fileID, '[PASSED] Dicretized Z-scores same as reference.\n');
    end
    
    % release some memory
    clear zScores;
    clear ref_zScores;
    
end


%% Check FC2doc for inference
% fprintf('\n--------[CHECK 2] Checking FC2doc results (factor composition inference):\n');
fprintf(log_fileID, '\n--------[CHECK 2] Checking FC2doc results (factor composition inference):\n');

% Check if dx1.dat and dx2.dat exist
if (~exist(fullfile(output_dir,'FC2doc','step1_output_inf_dx1.dat'), 'file'))
%     fprintf('[FAILED] Word counts document for ASD participants missing!\n');
    fprintf(log_fileID, '[FAILED] Word counts document for ASD participants missing!\n');
    error('[FAILED] Word counts document for ASD participants missing!');
end

if (~exist(fullfile(output_dir,'FC2doc','step1_output_dx2.dat'), 'file'))
%     fprintf('[FAILED] Word counts document for control participants missing!\n');
    fprintf(log_fileID, '[FAILED] Word counts document for control participants missing!\n');
    error('[FAILED] Word counts document for control participants missing!');
end

% Check discretized z scores, these are what were written into .dat word
% counts documents.
if (~exist(fullfile(output_dir,'FC2doc','step1_output_zScores.mat'), 'file'))
%     fprintf('[FAILED] Z scores missing!\n');
    fprintf(log_fileID, '[FAILED] Z scores missing!\n');
    error('[FAILED] Z scores missing!');
else
    zScores = load(fullfile(output_dir,'FC2doc','step1_output_inf_zScores.mat'));
    ref_zScores = load(fullfile(ref_dir,'FC2doc','step1_output_inf_zScores.mat'));
    
    if (~isequal(ref_zScores.discretized_z, zScores.discretized_z))
        diff_discretZ = max(max(abs(ref_zScores.discretized_z - zScores.discretized_z)));
%         fprintf('[FAILED] Dicretized Z-scores are different from reference, max abs diff = %f.\n]', diff_discretZ);
        fprintf(log_fileID, '[FAILED] Dicretized Z-scores differ from ref, max abs diff=%f.\n]', diff_discretZ);
        error('[FAILED] Dicretized Z-scores differ from ref, max abs diff=%f.]', diff_discretZ);
    else
%         fprintf('[PASSED] Dicretized Z-scores are the same as reference.\n');
        fprintf(log_fileID, '[PASSED] Dicretized Z-scores are the same as reference.\n');
    end
    
    % release some memory
    clear zScores;
    clear ref_zScores;
    
end

%% Check if final estimated model parameters are the same as reference
% fprintf('\n--------[CHECK 3] Checking if estimated model parameters are the same as reference:\n');
fprintf(log_fileID, '\n--------[CHECK 3] Checking if estimated model parameters same as reference:\n');

curr_output_dir = fullfile(output_dir,'estimate','k2','r1');
curr_ref_dir = fullfile(ref_dir,'estimate','k2','r1');

%%% beta estimate
if (~exist(fullfile(curr_output_dir,'final.beta'), 'file'))
%     fprintf('[FAILED] Final beta estimate result missing!\n');
    fprintf(log_fileID, '[FAILED] Final beta estimate result missing!\n');
    error('[FAILED] Final beta estimate result missing!');
else
    beta = load(fullfile(curr_output_dir,'final.beta'));
    ref_beta = load(fullfile(curr_ref_dir,'final.beta'));
    if (~isequal(beta, ref_beta))
        diff_beta = max(max(abs(beta - ref_beta)));
%         fprintf('[FAILED] Final beta estimate is different from reference, max abs diff: %f.\n', diff_beta);
        fprintf(log_fileID, '[FAILED] Final beta estimate differs from ref, max abs diff=%f.\n', diff_beta);
        error('[FAILED] Final beta estimate differs from ref, max abs diff=%f.', diff_beta);
    else
%         fprintf('[PASSED] Final beta estimate is the same as reference.\n');
        fprintf(log_fileID, '[PASSED] Final beta estimate is the same as reference.\n');
    end
    
    % release some memory
    clear beta;
    clear ref_beta;
    
end

%%% rho estimate
if (~exist(fullfile(curr_output_dir,'final.rho'), 'file'))
%     fprintf('[FAILED] Final rho estimate result missing!\n');
    fprintf(log_fileID, '[FAILED] Final rho estimate result missing!\n');
    error('[FAILED] Final rho estimate result missing!');
else
    rho = load(fullfile(curr_output_dir,'final.rho'));
    ref_rho = load(fullfile(curr_ref_dir,'final.rho'));
    if (~isequal(rho, ref_rho))
        diff_rho = max(max(abs(rho - ref_rho)));
%         fprintf('[FAILED] Final rho estimate is different from reference, max abs diff: %f.\n', diff_rho);
        fprintf(log_fileID, '[FAILED] Final rho estimate differs from ref, max abs diff=%f.\n', diff_rho);
        error('[FAILED] Final rho estimate differs from ref, max abs diff=%f.', diff_rho);
    else
%         fprintf('[PASSED] Final rho estimate is the same as reference.\n');
        fprintf(log_fileID, '[PASSED] Final rho estimate is the same as reference.\n');
    end
    
    % release some memory
    clear rho;
    clear ref_rho;
    
end
        
%%% gamma estimate
if (~exist(fullfile(curr_output_dir,'final.gamma'), 'file'))
%     fprintf('[FAILED] Final gamma estimate result missing!\n');
    fprintf(log_fileID, '[FAILED] Final gamma estimate result missing!\n');
    error('[FAILED] Final gamma estimate result missing!');
else
    gamma = load(fullfile(curr_output_dir,'final.gamma'));
    ref_gamma = load(fullfile(curr_ref_dir,'final.gamma'));
    if (~isequal(gamma, ref_gamma))
        diff_gamma = max(max(abs(gamma - ref_gamma)));
%         fprintf('[FAILED] Final gamma estimate is different from reference, max abs diff: %f.\n', diff_gamma);
        fprintf(log_fileID, '[FAILED] Final gamma estimate differs from ref, max abs diff=%f.\n', diff_gamma);
        error('[FAILED] Final gamma estimate differs from ref, max abs diff=%f.', diff_gamma);
    else
%         fprintf('[PASSED] Final gamma estimate is the same as reference.\n');
        fprintf(log_fileID, '[PASSED] Final gamma estimate is the same as reference.\n');
    end
    
    % release some memory
    clear gamma;
    clear ref_gamma;
    
end
        
%%% likelihood
if (~exist(fullfile(curr_output_dir,'likelihood.dat'), 'file'))
%     fprintf('[FAILED] Likelihood file missing!\n');
    fprintf(log_fileID, '[FAILED] Likelihood file missing!\n');
    error('[FAILED] Likelihood file missing!');
else
    likelihood = load(fullfile(curr_output_dir,'likelihood.dat'));
    ref_likelihood = load(fullfile(curr_ref_dir,'likelihood.dat'));
    if (~isequal(likelihood, ref_likelihood))
        diff_likelihood = max(max(likelihood - ref_likelihood));
%         fprintf('[FAILED] Likelihood is different from reference, max abs diff: %f.\n', diff_likelihood);
        fprintf(log_fileID, '[FAILED] Likelihood differs from ref, max abs diff=%f.\n', diff_likelihood);
        error('[FAILED] Likelihood differs from ref, max abs diff=%f.', diff_likelihood);
    else
%         fprintf('[PASSED] Likelihood is the same as reference.\n');
        fprintf(log_fileID, '[PASSED] Likelihood is the same as reference.\n');
    end
    
    % release some memory
    clear likelihood;
    clear ref_likelihood;
    
end

%% Check if E(RSFC patterns|Factor) & Pr(Factor|Participant) are the same as reference
fprintf(log_fileID, '\n--------[CHECK 4] Checking E(RSFC patterns|Factor) & Pr(Factor|Participant):\n');

curr_output_dir = fullfile(output_dir,'visualizeFactors','k2','r1');
curr_ref_dir = fullfile(ref_dir,'visualizeFactors','k2','r1');

k = 2;
%%% Check if E(RSFC pattern | Factor) is the same as reference,
%%% for each factor
for i = 1:k
    mean = load(fullfile(curr_output_dir, ['mean' num2str(i) '.mat']));
    ref_mean = load(fullfile(curr_ref_dir, ['mean' num2str(i) '.mat']));
    
    if (~isequal(mean.mean_corrmat, ref_mean.mean_corrmat))
        diff_mean = max(max(abs(mean.mean_corrmat - ref_mean.mean_corrmat)));
        fprintf(log_fileID, '[FAILED] RSFC patterns of factor %d differ from ref, max abs diff=%f.\n', i, diff_mean);
        error('[FAILED] RSFC patterns of factor %d differ from ref, max abs diff=%f.', i, diff_mean);
    else
        fprintf(log_fileID, '[PASSED] RSFC patterns of factor %d is the same as reference.\n', i);
    end
end
% release some memory
clear mean;
clear ref_mean;

%%% Check if Pr(Factor | Participant) is the same as reference
factorComp = load(fullfile(curr_output_dir,'factorComp.txt'));
ref_factorComp = load(fullfile(curr_ref_dir,'factorComp.txt'));
if (~isequal(factorComp, ref_factorComp))
    diff_factorComp = max(max(abs(factorComp - ref_factorComp)));

    fprintf(log_fileID, '[FAILED] Factor composition differs from ref, max abs diff=%f.\n', diff_factorComp);
    error('[FAILED] Factor composition differs from ref, max abs diff=%f.', diff_factorComp);
else
%     fprintf('[PASSED] Factor composition is the same as reference.\n');
    fprintf(log_fileID, '[PASSED] Factor composition is the same as reference.\n');
end
% release some memory
clear factorComp;
clear ref_factorComp;

%% Check if inferred factor composition is the same as reference
% fprintf('\n--------[CHECK 5] Checking if inferred factor composition is the same:\n');
fprintf(log_fileID, '\n--------[CHECK 5] Checking if inferred factor composition is the same:\n');

curr_output_dir = fullfile(output_dir, 'inference');
curr_ref_dir = fullfile(ref_dir, 'inference');

%%% Factor composition
factorComp_files = dir(fullfile(curr_output_dir, '*k2r1_factorComp.txt'));
if isempty(factorComp_files)
%     fprintf('[FAILED] Factor composition text file missing!\n');
    fprintf(log_fileID, '[FAILED] Factor composition text file missing!\n');
    error('[FAILED] Factor composition text file missing!');
else
    factorComp = load(fullfile(curr_output_dir, factorComp_files(1).name));
    ref_factorComp_files = dir(fullfile(curr_ref_dir, '*k2r1_factorComp.txt')) ;
    ref_factorComp = load(fullfile(curr_ref_dir, ref_factorComp_files(1).name));
    
    if (~isequal(factorComp, ref_factorComp))
        diff_factorComp = max(max(abs(factorComp - ref_factorComp)));
%         fprintf('[FAILED] Factor composition is different 
        fprintf(log_fileID, '[FAILED] Factor composition differs from ref, max abs diff=%f.\n', diff_factorComp);
        error('[FAILED] Factor composition differs from ref, max abs diff=%f.', diff_factorComp);
    else
%         fprintf('[PASSED] Factor composition is the same as the reference.\n');
        fprintf(log_fileID, '[PASSED] Factor composition is the same as the reference.\n');
    end
    
    % release some memory
    clear factorComp;
    clear ref_factorComp;
    
end

%% Close log file
fprintf(log_fileID,'\n\n');
fclose(log_fileID);

end

