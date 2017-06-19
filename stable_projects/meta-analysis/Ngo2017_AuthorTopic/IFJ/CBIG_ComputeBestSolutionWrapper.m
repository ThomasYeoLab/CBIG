function CBIG_ComputeBestSolutionWrapper(allKs, seeds)
  ALPHA     = 100;
  ETA       = 0.01;
  
  for K = allKs
    estimatesDir = 'outputs';
    CBIG_ComputeBestSolution(estimatesDir, K, seeds, ALPHA, ETA);
  end