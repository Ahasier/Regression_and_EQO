function computeAndSaveRegressionResults(numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod, meshGrid, varargin)
% COMPUTEANDSAVEREGRESSIONRESULTS Computes and saves the regression results.
% 
% INPUTS:
%   numPermutations: Number of random permutations.
%   phylogenyDependency: 0 if beta's are independent of phylogeny.
%   noiseLevel: Noise level for regression.
%   numGroups: Number of taxa in a functional group.
%   numSamples: Number of samples in the trainning data.
%   regressionMethod: Regression method in use.
%   varargin: Option pairs. i.e.
%   'BetaEps',0.5,'Beta0',1,'Threshold','cv','RealAbd','On'.
%
% Example usage: computeAndSaveRegressionResults(10, 0, 1, 5, 20, 'LASSO', 'maxLambda', 10, 'Threshold', 'cv');

% Add necessary paths
initializePaths();

% Extracting settings and paths
[treeData, settings, fullIdentifier, resultsPath, betaResultsPath] = setupPathsAndSettings(varargin{:});

% Handle based on real or synthetic data
if usingRealData(settings)
    results = processRealData(settings, treeData, numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod);
else
    results = processSyntheticData(settings, treeData, numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod);
end

% Save results to disk
saveResults(results, resultsPath, betaResultsPath, regressionMethod, fullIdentifier, numberOfTaxaInAGroup, numSamples, meshGrid);
end

%% Helper functions
function initializePaths()
functionsPath = [pwd, '/functions'];
addpath(functionsPath);
end

function [treeData, settings, fullIdentifier, resultsPath, betaResultsPath] = setupPathsAndSettings(varargin)
resultsPath = 'results/';
betaResultsPath = 'results/Betas/';
load('tree100taxa.mat');
treeData = tree100;
[settings, fullIdentifier] = setOptionsAndNames(varargin{:});
end

function results = processRealData(settings, treeData, numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod)
% Generate Data
[abundanceData, functionalOutput] = generateData(treeData, phylogenyDependency, numberOfTaxaInAGroup, noiseLevel, numSamples, settings);

% Compute regression and cross-validation over different thresholds on data
[allCoefficients, allErrorsAtAllThresholds, crossValidatedCoefficients, crossValidateThreshold, coefficientsStdDev, OutSampleError] = computeRegressionAndCrossValidation(abundanceData, functionalOutput, numPermutations, regressionMethod, settings.Beta0);

% Return results as a structure
results = struct('abundanceData', abundanceData, 'functionalOutput', functionalOutput, 'crossValidatedCoefficients', crossValidatedCoefficients, 'allCoefficients', allCoefficients, 'allErrorsAtAllThresholds', allErrorsAtAllThresholds, 'crossValidateThreshold', crossValidateThreshold, 'coefficientsStdDev', coefficientsStdDev, 'MeanSquaredErrorOutOfSample', OutSampleError);
end

function results = processSyntheticData(settings, treeData, numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod)
% Generate Data
[abundanceData, functionalOutput, syntheticCoefficients] = generateData(treeData, phylogenyDependency, numberOfTaxaInAGroup, noiseLevel, numSamples, settings);

% Compute regression and cross-validation over different thresholds on data
[allCoefficients, allErrorsAtAllThresholds, crossValidatedCoefficients, crossValidateThreshold, coefficientsStdDev, OutSampleError] = computeRegressionAndCrossValidation(abundanceData, functionalOutput, numPermutations, regressionMethod, settings.Beta0);

% Calculate Accuracy
accuracy = calculateAccuracy(crossValidatedCoefficients, syntheticCoefficients, settings);

% Return results as a structure
results = struct('abundanceData', abundanceData, 'functionalOutput', functionalOutput, 'syntheticCoefficients', syntheticCoefficients, 'crossValidatedCoefficients', crossValidatedCoefficients, 'allCoefficients', allCoefficients, 'allErrorsAtAllThresholds', allErrorsAtAllThresholds, 'crossValidateThreshold', crossValidateThreshold, 'coefficientsStdDev', coefficientsStdDev, 'MeanSquaredErrorOutOfSample', OutSampleError,'accuracy', accuracy);
end