function [numPermutations, phylogenyDependency, noiseLevel, meshGrid, settings] = initializations(regressionMethod)
% Initialize paths
initializePaths();

% Set parameters from JSON file
paramsFilename = 'configurations/basicParams.json';
[numPermutations, phylogenyDependency, noiseLevel, meshGrid] = setParams(paramsFilename);

% Set other parameters using setOptionsAndNames function
settings = setOptionsAndNames();

% Define global variable paths for where to load data or store results
global paths
paths = SetPathsForDataAndResults('data', 'results', 'betaResults','accuracyResults', 'tcmResults');

% If using EQO, add the path to `Rscript` to PATH in MATLAB
if strcmp(regressionMethod, 'EQO')
    setenv('PATH', [getenv('PATH') ':/usr/local/bin/']);
end
end