function CBIG_factorHierachy(k_bestRun)

% CBIG_factorHierachy(k_bestRun)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

noSplits = size(k_bestRun, 1)-1;

K_0 = k_bestRun{1, 1};
K = k_bestRun{end, 1}-1;

which2_order_avgCorr_beta = cell(noSplits, 3);
which2_order_avgCorr_gamma = cell(noSplits, 3);

% k-to-(k+1) split
for k = K_0:K
    dir = k_bestRun{cell2mat(k_bestRun(:, 1))==k, 2};
    kPlus1 = k+1;
    dirPlus1 = k_bestRun{cell2mat(k_bestRun(:, 1))==kPlus1, 2};
    
    %------------------------- Beta
    kBeta = load([dir 'final.beta']);
    kBeta = exp(kBeta);
    kPlus1Beta = load([dirPlus1 'final.beta']);
    kPlus1Beta = exp(kPlus1Beta);
    
    %------------------------- Gamma
    kGamma = load([dir 'final.gamma']);
    kGamma = bsxfun(@times, kGamma, 1./(sum(kGamma, 2))); % normalization
    kPlus1Gamma = load([dirPlus1 'final.gamma']);
    kPlus1Gamma = bsxfun(@times, kPlus1Gamma, 1./(sum(kPlus1Gamma, 2))); % normalization
    
    % kPlus1-choose-2 possible merges
    toCombine = combnk(1:kPlus1, 2);
    
    % Exhaustive search of pairs to combine
    noComb = size(toCombine, 1);
    newKbeta_orders_costs = cell(noComb, 3);
    newKgamma_orders_costs = cell(noComb, 3);
    
    % For every possible merge
    for idx = 1:noComb
        idx1 = toCombine(idx, 1);
        idx2 = toCombine(idx, 2);
        
        % Derive another k-subtype estimate from the (k+1)-subtype one
        mean_row = (kPlus1Beta(idx1, :)+kPlus1Beta(idx2, :))/2;
        newKbeta = [ ...
            mean_row; ...
            kPlus1Beta(1:(idx1-1), :); ...
            kPlus1Beta((idx1+1):(idx2-1), :); ...
            kPlus1Beta((idx2+1):end, :)];
        newKbeta_orders_costs{idx, 1} = newKbeta;
        sum_col = kPlus1Gamma(:, idx1)+kPlus1Gamma(:, idx2);
        newKgamma = [ ...
            sum_col ...
            kPlus1Gamma(:, 1:(idx1-1)) ...
            kPlus1Gamma(:, (idx1+1):(idx2-1)) ...
            kPlus1Gamma(:, (idx2+1):end)];
        newKgamma_orders_costs{idx, 1} = newKgamma;
        
        %%% Reorder subtypes to obtain the maximal correlation coefficients
        % Construct the COST matrix
        %         pos1 pos2 ...
        % topic1
        % topic2
        % ...
        costMatBeta = zeros(k, k);
        costMatGamma = zeros(k, k);
        for rowIdx = 1:k
            for colIdx = 1:k
                % Assign kBeta (jobs, column) to newKbeta (workers, row)
                corrMatBeta = corrcoef(newKbeta(rowIdx, :)', kBeta(colIdx, :)');
                costMatBeta(rowIdx, colIdx) = 1-corrMatBeta(1, 2);
                corrMatGamma = corrcoef(newKgamma(:, rowIdx), kGamma(:, colIdx));
                costMatGamma(rowIdx, colIdx) = 1-corrMatGamma(1, 2);
            end
        end
        
        % Run the Hungarian matching algorithm
        % order: each row's matched column
        [orderBeta, costBeta] = munkres(costMatBeta);
        newKbeta_orders_costs{idx, 2} = orderBeta;
        newKbeta_orders_costs{idx, 3} = costBeta;
        [orderGamma, costGamma] = munkres(costMatGamma);
        newKgamma_orders_costs{idx, 2} = orderGamma;
        newKgamma_orders_costs{idx, 3} = costGamma;
    end
    
    % Find the optimal pairs to combine and the order
    [~, bestIdxBeta] = min(cell2mat(newKbeta_orders_costs(:, 3)));
    [~, bestIdxGamma] = min(cell2mat(newKgamma_orders_costs(:, 3)));
    
     % Check if same split but different ordering are found for two probabilities
    if bestIdxBeta ~= bestIdxGamma
        fprintf('Different merges and orderings are found for the %i-to-%i split\n', k, k+1);

        str = input('Accomodate beta (b) or gamma (g) or neither (n)? b/g/n: ', 's');
        if str == 'b'
            bestIdxGamma = bestIdxBeta;
        elseif str == 'g'
            bestIdxBeta = bestIdxGamma;
        elseif str == 'n'
            % Keep as is
        else
            error('Wrong input!');
        end
    end
    
    bestNewKbeta = newKbeta_orders_costs{bestIdxBeta, 1};
    bestOrderBeta = newKbeta_orders_costs{bestIdxBeta, 2};
    bestNewKgamma = newKgamma_orders_costs{bestIdxGamma, 1};
    bestOrderGamma = newKgamma_orders_costs{bestIdxGamma, 2};
    
    % Reorder
    reorderdedBestNewKbeta = bestNewKbeta;
    reorderdedBestNewKgamma = bestNewKgamma;
    for idx = 1:size(reorderdedBestNewKbeta, 1)
        reorderdedBestNewKbeta(bestOrderBeta(idx), :) = bestNewKbeta(idx, :);
        reorderdedBestNewKgamma(:, bestOrderGamma(idx)) = bestNewKgamma(:, idx);
    end
    
    corrSumBeta = 0;
    corrSumGamma = 0;
    for idx = 1:k
        corrMatBeta = corrcoef(reorderdedBestNewKbeta(idx, :)', kBeta(idx, :)');
        corrMatGamma = corrcoef(reorderdedBestNewKgamma(:, idx), kGamma(:, idx));
        corrSumBeta = corrSumBeta+corrMatBeta(1, 2);
        corrSumGamma = corrSumGamma+corrMatGamma(1, 2);
    end
    
    which2_order_avgCorr_beta{k-1, 1} = toCombine(bestIdxBeta, :);
    which2_order_avgCorr_beta{k-1, 2} = bestOrderBeta;
    which2_order_avgCorr_beta{k-1, 3} = corrSumBeta/k;
    
    which2_order_avgCorr_gamma{k-1, 1} = toCombine(bestIdxGamma, :);
    which2_order_avgCorr_gamma{k-1, 2} = bestOrderGamma;
    which2_order_avgCorr_gamma{k-1, 3} = corrSumGamma/k;
end

%------------------------ Plot

CBIG_plotSetup;

close all;
figure('Position', [100, 600, 400, 800]);

xTickLabels = {};
for idx = K_0:K
    xTickLabels = [xTickLabels, [num2str(idx) '-' num2str(idx+1)]];
end

for idx = 1:2
    subplot(2, 1, idx)
    if idx == 1
        yStr = 'Average Correlation';
        titleStr = '(A) Pr(Voxel | Factor)'
        avgCorr = cell2mat(which2_order_avgCorr_beta(:, 3))
    else
        yStr = 'Average Correlation';
        titleStr = '(B) Pr(Factor | Patient)'
        avgCorr = cell2mat(which2_order_avgCorr_gamma(:, 3))
    end
    plot(K_0:K, avgCorr, 'bo-', 'LineWidth', 3, 'MarkerSize', 10);
    grid on;
    ax = gca;
    set(ax, 'XLim', [K_0, K]);
    set(ax, 'XTick', K_0:K);
    set(ax, 'XTickLabel', xTickLabels);
    set(ax, 'YLim', [0, 1]);
    set(ax, 'YTick', 0:0.1:1);
    set(ax, 'FontSize', 20);
    xlabel('#Factors');
    ylabel(yStr);
    title(titleStr, 'fontsize', 20, 'fontweight', 'bold');
    box off;
end