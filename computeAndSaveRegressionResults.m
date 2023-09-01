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

% If on diagnostic mod, save results for diagnostics as well; otherwise
% only save the main output results.
if onDiagnosticMod(settings)
    % Handle based on real or synthetic data
    if usingRealData(settings)
        [results, resultsForDiagnostics] = processRealData(settings, treeData, numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod);
    else
        [results, resultsForDiagnostics] = processSyntheticData(settings, treeData, numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod);
    end
    
    % Combine main output results with those for diagnostic uses.
    results = combineResults(results, resultsForDiagnostics);
else
    % Handle based on real or synthetic data
    if usingRealData(settings)
        results = processRealData(settings, treeData, numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod);
    else
        results = processSyntheticData(settings, treeData, numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod);
    end
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

function [results, resultsForDiagnostics] = processRealData(settings, treeData, numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod)
% Generate Data
[abundanceData, functionalOutput] = generateData(treeData, phylogenyDependency, numberOfTaxaInAGroup, noiseLevel, numSamples, settings);

% Compute regression and cross-validation over different thresholds on data
[allCoefficients, allBestCoefficients, crossValidatedCoefficients, allOutSampleErrorsAtAllThresholds, crossValidatedOutSampleError, allOutSampleR2AtAllThresholds, allOptimalThresholds, coefficientsStdDev] = computeRegressionAndCrossValidation(abundanceData, functionalOutput, numPermutations, regressionMethod, settings.Beta0);

% Return results as a structure
results = struct('abundanceData', abundanceData, 'functionalOutput', functionalOutput, 'crossValidatedCoefficients', crossValidatedCoefficients, 'coefficientsStdDev', coefficientsStdDev, 'MeanSquaredErrorOutOfSample', crossValidatedOutSampleError);
resultsForDiagnostics = struct('allCoefficients', allCoefficients, 'allBestCoefficients', allBestCoefficients, 'allErrorsAtAllThresholds', allOutSampleErrorsAtAllThresholds, 'allOutSampleR2AtAllThresholds', allOutSampleR2AtAllThresholds, 'allOptimalThresholds', allOptimalThresholds);
end

function [results, resultsForDiagnostics] = processSyntheticData(settings, treeData, numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod)
% Generate Data
[abundanceData, functionalOutput, syntheticCoefficients] = generateData(treeData, phylogenyDependency, numberOfTaxaInAGroup, noiseLevel, numSamples, settings);

% Compute regression and cross-validation over different thresholds on data
[allCoefficients, allBestCoefficients, crossValidatedCoefficients, allOutSampleErrorsAtAllThresholds, crossValidatedOutSampleError, allOutSampleR2AtAllThresholds, allOptimalThresholds, coefficientsStdDev] = computeRegressionAndCrossValidation(abundanceData, functionalOutput, numPermutations, regressionMethod, settings.Beta0);

% Calculate Accuracy
accuracy = calculateAccuracy(crossValidatedCoefficients, syntheticCoefficients, settings);

% Return results as a structure
results = struct('abundanceData', abundanceData, 'functionalOutput', functionalOutput, 'syntheticCoefficients', syntheticCoefficients, 'crossValidatedCoefficients', crossValidatedCoefficients, 'coefficientsStdDev', coefficientsStdDev, 'MeanSquaredErrorOutOfSample', crossValidatedOutSampleError, 'accuracy', accuracy);
resultsForDiagnostics = struct('allCoefficients', allCoefficients, 'allBestCoefficients', allBestCoefficients, 'allErrorsAtAllThresholds', allOutSampleErrorsAtAllThresholds, 'allOutSampleR2AtAllThresholds', allOutSampleR2AtAllThresholds, 'allOptimalThresholds', allOptimalThresholds);
end

function isDiagnosticMod = onDiagnosticMod(settings)
isDiagnosticMod = isfield(settings, 'DiagnosticMod') && strcmp(settings.DiagnosticMod, 'On');
end

function sCombined = combineResults(s1, s2)
sCombined = s1;
fields = fieldnames(s2);
for i = 1:numel(fields)
    sCombined.(fields{i}) = s2.(fields{i});
end
end