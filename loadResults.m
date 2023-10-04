function results = loadResults(regressionMethod, fullIdentifier, numberOfTaxaInAGroup, numSamples)
% Recall paths for results storage
global paths

% Construct the file path where other simulation results will be saved
betaResultsFilePath = [paths.betaResultsPath, 'Betas_', regressionMethod, fullIdentifier, '_K', num2str(numberOfTaxaInAGroup), '_nSpl', num2str(numSamples), '.mat'];

% Save various results to a MAT file for later analysis
load(betaResultsFilePath, 'results');
end