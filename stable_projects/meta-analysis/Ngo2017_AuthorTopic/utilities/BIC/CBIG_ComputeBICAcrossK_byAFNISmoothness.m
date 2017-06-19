function CBIG_ComputeBICAcrossK_byAFNISmoothness(allKs, solution_dir, smoothness_dir, solution_type, output_dir)
  allBICs = zeros(length(allKs), 1);
  allLogLikelihoods = zeros(length(allKs), 1);
  allComplexityCosts = zeros(length(allKs), 1);
  allFreeDimensionsCounts = zeros(length(allKs), 1);
  count = 1;
  
  bic_dir = fullfile(output_dir, solution_type, 'BIC_3dFWHMx');
  system(['mkdir -p ' bic_dir]);
  
  for K = allKs  
    disp(['K: ' num2str(K)]);
    smoothness = load(fullfile(smoothness_dir, ['K' num2str(K) '_smoothness.mat']));
    params_data = load(fullfile(solution_dir, [solution_type '_K' num2str(K, '%03d') '.mat']));
    
    freeDimensionsCount = smoothness.sum_afni_resels_count + size(params_data.params.theta, 1) * (K-1);
    [bic, logLikelihood, complexityCost] = CBIG_ComputeBICOnData(params_data.params, freeDimensionsCount);
    
    save(fullfile(bic_dir, ['K' num2str(K) '_BIC.mat']), 'bic');
    allBICs(count) = bic;
    allLogLikelihoods(count) = logLikelihood;
    allFreeDimensionsCounts(count) = freeDimensionsCount;
    allComplexityCosts(count) = complexityCost;
    count = count + 1;
  end

  figure;
  plot(allKs, allBICs, '-ko', 'LineWidth', 2);
  title('Bayesian Information Criteria (BIC)', 'FontSize', 19);
  xlabel('# Components', 'FontSize', 19);
  ylabel('Bayesian Information Criterion (BIC)', 'FontSize', 19);
  set(gca, 'XTick', allKs);
  set(gca, 'FontSize', 19);
  set(gca, 'TickDir', 'out');
  box off;
  num_ticks = 4;
  L = get(gca, 'YLim');
  set(gca, 'YTick', linspace(L(1), L(2), num_ticks));
  saveas(gcf, fullfile(bic_dir, 'BIC_3dFWHMx'), 'epsc');
