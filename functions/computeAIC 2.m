function aicValue = computeAIC(n, numTaxa, trainingData, trainingOutput, sortedTaxaIndices)
idx = sortedTaxaIndices(n + 1:end);
groupAssemblage = ones(numTaxa, 1);
%     groupAssemblage = coefficients;
groupAssemblage(idx) = 0;

logLikelihood = regressionLikelihood(trainingData, trainingOutput, groupAssemblage);

%     [r2, ssr] = computeRSquared(trainingData, trainingOutput, groupAssemblage);
%     L = exp(-ssr/2);
aic = 2*n - 2*logLikelihood;

% Store the AIC value
aicValue = aic;
end