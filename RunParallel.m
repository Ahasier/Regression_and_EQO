function results = RunParallel(numberOfTaxaInAGroup, numSamples, regressionMethod)
% Initialize all neccessary parametters from configurations files
[numPermutations, phylogenyDependency, noiseLevel, meshGrid, settings, fullIdentifier] = initializations(regressionMethod);

% Run computeAndSaveRegressionResults for the given numberOfTaxaInAGroup and numSamples
results = computeAndSaveRegressionResults(numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod, meshGrid, settings, fullIdentifier);
end
