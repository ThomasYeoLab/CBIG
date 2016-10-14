function T_S_C_p = CBIG_blAge_edu_amyloid_apoe_disDura_onsetAge(subInfoFile, GMICVFile, K, Q, dx)

% T_S_C_p = CBIG_blAge_edu_amyloid_apoe_disDura_onsetAge(subInfoFile, GMICVFile, K, Q, dx)
% Which quantity? Q
% 1 for baseline age, 2 for edu, 3 for amyloid, 4 for APOE-2, 5 for APOE-4,
% 6 for disease duration, 7 for onset age
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


T_S_C_p = cell(1, 4);

%% Baseline info

[rid_dx, ~, rid_prob, rid_age, ~, rid_edu] = get_data(subInfoFile, GMICVFile, K);

% RIDs of 188 AD, 147 a+ MCI, 43 a+ AD
if dx == 3 % 188 AD
    rid = rid_dx(rid_dx(:, 2)==3, 1);
elseif dx == 2 % a+ MCI
    rid = CBIG_select_predementedGrp(rid_dx, 'a+_mci_147');
elseif dx == 1 % a+ CN
    rid = CBIG_select_predementedGrp(rid_dx, 'a+_hc_43');
end

%%

if Q == 1 % baseline age
    ind = ismember(rid_age(:, 1), rid);
    rid_y = rid_age(ind, :);
    % Factor probability
    ind = ismember(rid_prob(:, 1), rid_y(:, 1));
    rid_prob = rid_prob(ind, :);
    % Check if RID matched
    assert(isequal(rid_y(:, 1), rid_prob(:, 1)));
    qName = '************* BASELINE AGE *************';
elseif Q == 2 % edu
    ind = ismember(rid_edu(:, 1), rid);
    rid_y = rid_edu(ind, :);
    % Factor probability
    ind = ismember(rid_prob(:, 1), rid_y(:, 1));
    rid_prob = rid_prob(ind, :);
    % Check if RID matched
    assert(isequal(rid_y(:, 1), rid_prob(:, 1)));
    qName = '************* EDU *************';
elseif Q == 3 % amyloid
    amyloid = CBIG_get_amyloid(rid);
    rid_y = cell2mat(amyloid(:, [1 4]));
    % Factor probability: not everyone has amyloid results
    ind = ismember(rid_prob(:, 1), rid_y(:, 1));
    rid_prob = rid_prob(ind, :);
    % Check if RID matched
    assert(isequal(rid_y(:, 1), rid_prob(:, 1)));
    qName = '************* AMYLOID *************';
elseif Q == 4 % ApoE-2
    apoe = CBIG_get_apoe(rid);
    rid_y = cell2mat(apoe(:, [1 3]));
    % Factor probability
    ind = ismember(rid_prob(:, 1), rid_y(:, 1));
    rid_prob = rid_prob(ind, :);
    % Check if RID matched
    assert(isequal(rid_y(:, 1), rid_prob(:, 1)));
    qName = '************* APOE 2 *************';
elseif Q == 5 % ApoE-4
    apoe = CBIG_get_apoe(rid);
    rid_y = cell2mat(apoe(:, [1 4]));
    % Factor probability
    ind = ismember(rid_prob(:, 1), rid_y(:, 1));
    rid_prob = rid_prob(ind, :);
    % Check if RID matched
    assert(isequal(rid_y(:, 1), rid_prob(:, 1)));
    qName = '************* APOE 4 *************';
elseif Q == 6 % AD only: disease duration
    rid_durations = CBIG_compute_disDura(subInfoFile, GMICVFile, K);
    ind = ismember(rid_durations(:, 1), rid);
    rid_y = rid_durations(ind, :);
    % Factor probability
    ind = ismember(rid_prob(:, 1), rid_y(:, 1));
    rid_prob = rid_prob(ind, :);
    % Check if RID matched
    assert(isequal(rid_y(:, 1), rid_prob(:, 1)));
    qName = '************* YEARS SINCE ONSET *************';
elseif Q == 7 % AD only: onset age
    rid_durations = CBIG_compute_disDura(subInfoFile, GMICVFile, K);
    % disDura = yyyyBl - yyyyOnset
    % ageBl = yyyyyBl - yyyyBirth
    % ageOnset = yyyyOnset - yyyyBirth
    %          = ageBl - disDura
    ind = ismember(rid_age(:, 1), rid_durations(:, 1));
    rid_ageBl = rid_age(ind, :);
    rid_y = [rid_ageBl(:, 1) rid_ageBl(:, 2)-rid_durations(:, 2)];
    % Factor probability
    ind = ismember(rid_prob(:, 1), rid_y(:, 1));
    rid_prob = rid_prob(ind, :);
    % Check if RID matched
    assert(isequal(rid_y(:, 1), rid_prob(:, 1)));
    qName = '************* ONSET AGE *************';
else
    error('Wrong selection')
end

y = rid_y(:, 2);

disp(qName);
fprintf('\n');
fprintf('# of responses: %d\n', numel(y));

%% Compute stats

[wMean, wStdDev] = CBIG_compute_weightedMeanStdDev(rid_prob(:, 2:end), y);

if K == 3
    % What is the subtype order?
    % Offset by 1, the RID column
    TEM_COL = 1+2;
    SUB_COL = 1+3;
    COR_COL = 1+1;
    fprintf('\n');
    fprintf('!!! Subtype order assumed to be: cortical, temporal, subcortical.\n');
    fprintf('!!! Make sure this is true in the gamma file.\n');
    %------------------------------- Hypothesis testing
    X = rid_prob(:, [TEM_COL SUB_COL]);
    [b, ~, stats] = glmfit(X, y);
    %---------------- Overall
    fprintf('\n');
    disp('Overall');
    % For linear regression, the p-values are exact.
    % H*b = c
    H = [0 1 -1; 0 1 0];
    c = [0; 0];
    p = linhyptest(b, stats.covb, c, H, stats.dfe)
    %---------------- Pairwise comparisons
    % T - S
    H = [0 1 -1];
    c = 0;
    if H*b > 0
        fprintf('T > S\n');
    else
        fprintf('T < S\n');
    end
    linhyptest(b, stats.covb, c, H, stats.dfe)
    % T - C
    H = [0 1 0];
    c = 0;
    if H*b > 0
        fprintf('T > C\n');
    else
        fprintf('T < C\n');
    end
    linhyptest(b, stats.covb, c, H, stats.dfe)
    % S - C
    H = [0 0 1];
    c = 0;
    if H*b > 0
        fprintf('S > C\n');
    else
        fprintf('S < C\n');
    end
    linhyptest(b, stats.covb, c, H, stats.dfe)
else
    error('Unconfigured; refer to older versions.');
end

T_S_C_p{1} = sprintf('%.3f (%.3f)', wMean(TEM_COL-1), wStdDev(TEM_COL-1));
T_S_C_p{2} = sprintf('%.3f (%.3f)', wMean(SUB_COL-1), wStdDev(SUB_COL-1));
T_S_C_p{3} = sprintf('%.3f (%.3f)', wMean(COR_COL-1), wStdDev(COR_COL-1));
T_S_C_p{4} = sprintf('%e', p);