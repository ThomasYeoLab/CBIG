function [bic, logLikelihood, complexityCost] = ComputeBICOnData(params, freeDimensionsCount)
  logLikelihood = 0;
  wordCountSum = 0;
  
  beta = params.beta;
  theta = params.theta;

  for d = 1:params.D
    wordIndices = logical(params.Wb(d, :));
    betaForDoc = beta(:, wordIndices);
    wordCount = params.Wc{d}';

    authorsForDoc = logical(full(params.doc_by_author(d, :)));
    authorsForDocCount = sum(authorsForDoc);
    thetaForDoc = theta(authorsForDoc, :);
    
    docLogLikelihood = betaForDoc' * thetaForDoc';
    docLogLikelihood = log(sum(docLogLikelihood, 2));
    docLogLikelihood = bsxfun(@times, wordCount', docLogLikelihood);
    docLogLikelihood = sum(docLogLikelihood) - log(authorsForDocCount) * sum(wordCount);
    
    logLikelihood = logLikelihood + docLogLikelihood;
    wordCountSum = wordCountSum + sum(wordCount);
  end

  bic = logLikelihood - 0.5 * freeDimensionsCount * log(wordCountSum);
  complexityCost = -0.5 * freeDimensionsCount * log(wordCountSum);