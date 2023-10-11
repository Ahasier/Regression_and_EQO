clc; clear;

% Initialize paths
initializePaths();

% Set parameters from JSON file
paramsFilename = 'configurations/basicParams.json';
[numPermutations, phylogenyDependency, noiseLevel, meshGrid] = setParams(paramsFilename);

regressionMethod = 'EQO';

% Set other parameters using setOptionsAndNames function
[settings, fullIdentifier] = setOptionsAndNames();

% Define global variable paths for where to load data or store results
global paths
paths = SetPathsForDataAndResults('data', 'results', 'betaResults','accuracyResults', 'tcmResults');

numberOfTaxaInAGroup = 10;
numSamples = 80;

results = loadResults(regressionMethod, fullIdentifier, numberOfTaxaInAGroup, numSamples);

maxImportanceValue = max(results.importanceValues)
compareCoefficients = [results.importanceValues, results.crossValidatedCoefficients, results.syntheticCoefficients]
