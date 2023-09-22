function aicValue = computeAIC(n, trainingData, trainingOutput, coefficients, sortedTaxaIndices, settings)
idx = sortedTaxaIndices(n + 1:end);
% groupAssemblage = ones(numTaxa, 1);
groupAssemblage = coefficients;
groupAssemblage(idx) = 0;

if isfield(settings, 'weight') && strcmp(settings.weight, 'On')
    % Loop over different weights to pick an optimal one
    weights = setWeightsRange();
    lenWeights = length(weights);
    aics = zeros(1, lenWeights);
    
    for w = 1:lenWeights
        weight = weights(w);
        groupAssemblage = weight.*groupAssemblage;
        logLikelihood = regressionLikelihood(trainingData, trainingOutput, groupAssemblage);
        aics(w) = 2*n - 2*logLikelihood;
    end
    
    aicValue = min(aics);
else
    logLikelihood = regressionLikelihood(trainingData, trainingOutput, groupAssemblage);
    
    aicValue = 2*n - 2*logLikelihood;
end
end