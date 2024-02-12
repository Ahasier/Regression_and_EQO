function accuracies = runWithVaryingEpsilon(numberOfTaxaInAGroup, numSamples, regressionMethod, betaEps_list)
% Initialize all neccessary parametters from configurations files
[numPermutations, phylogenyDependency, noiseLevel, meshGrid, settings] = initializations(regressionMethod);

% Create a full identifier string based on the options
fullIdentifier = createIdentifier(settings);

% Run computeAndSaveRegressionResults with different BetaEps
for n = 1:length(betaEps_list)
    settings.BetaEps = betaEps_list(n);
    % Run computeAndSaveRegressionResults for the given numberOfTaxaInAGroup and numSamples
    results = computeAndSaveRegressionResults(numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod, meshGrid, settings, fullIdentifier);
    accuracies(n) = results.accuracy;
end

% Recall paths for results storage
global paths
% Construct the file path where other simulation results will be saved
accResultsFilePath = [paths.accuracyResultsPath, 'accuraciesVaryingEps_', regressionMethod, fullIdentifier, '_K', num2str(numberOfTaxaInAGroup), '_nSpl', num2str(numSamples), '.mat'];
save(accResultsFilePath, 'accuracies');
end