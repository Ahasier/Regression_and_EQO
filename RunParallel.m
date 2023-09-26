function results = RunParallel(numberOfTaxaInAGroup, numSamples, regressionMethod)
% Initialize paths
initializePaths();

% Set parameters from JSON file
paramsFilename = 'configurations/basicParams.json';
[numPermutations, phylogenyDependency, noiseLevel, meshGrid] = setParams(paramsFilename);

% Set other parameters using setOptionsAndNames function
[settings, fullIdentifier] = setOptionsAndNames();

% Define global variable paths for where to load data or store results
global paths
paths = SetPathsForDataAndResults('data', 'results', 'betaResults','accuracyResults', 'tcmResults');

% If using EQO, add the path to `Rscript` to PATH in MATLAB
if strcmp(regressionMethod, 'EQO')
    setenv('PATH', [getenv('PATH') ':/usr/local/bin/']);
end

% Run computeAndSaveRegressionResults for the given numberOfTaxaInAGroup and numSamples
results = computeAndSaveRegressionResults(numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod, meshGrid, settings, fullIdentifier);
end
