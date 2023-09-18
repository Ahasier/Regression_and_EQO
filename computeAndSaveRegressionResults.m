function results = computeAndSaveRegressionResults(numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod, meshGrid, settings, fullIdentifier)
% COMPUTEANDSAVEREGRESSIONRESULTS Computes and saves the regression results.
% 
% INPUTS:
%   numPermutations: Number of random permutations for the analysis.
%   phylogenyDependency: Dependency parameter related to phylogeny (0 if beta's are independent of phylogeny).
%   noiseLevel: Noise level to consider for the regression.
%   numberOfTaxaInAGroup: Number of taxa grouped together in a functional group.
%   numSamples: Number of samples in the training data.
%   regressionMethod: Chosen regression method (e.g., 'LASSO').
%
% OUTPUTS:
%   results: Regression results based on the specified parameters.
%
% Example usage: computeAndSaveRegressionResults(10, 0, 1, 5, 20, 'LASSO', 'maxLambda', 10, 'Threshold', 'cv');

% Load tree data
treeData = loadTreeData();

% Process data to get results
results = processData(settings, treeData, numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod);

% Save results to disk
saveResults(results, regressionMethod, fullIdentifier, numberOfTaxaInAGroup, numSamples, meshGrid);
end

%% Helper functions
% Load tree data
function treeData = loadTreeData()
global paths

load([paths.data,'tree100taxa.mat'], 'tree100');
treeData = tree100;
end

function results = processData(settings, treeData, numPermutations, phylogenyDependency, noiseLevel, numberOfTaxaInAGroup, numSamples, regressionMethod)
% 1. Generate Data
[abundanceData, functionalOutput, syntheticCoefficients] = generateData(treeData, phylogenyDependency, numberOfTaxaInAGroup, noiseLevel, numSamples, settings);

% 2. Compute regression and cross-validation over different thresholds on data
[TCM, importanceValues, crossValidatedCoefficients] = computeRegressionAndCrossValidation(abundanceData, functionalOutput, numPermutations, regressionMethod, settings);

% 3. If usimg mock data, calculate accuracy
if ~usingRealData(settings)
    accuracy = calculateAccuracy(crossValidatedCoefficients, syntheticCoefficients, settings);
end

% 4. Return results as a structure
results = struct('abundanceData', abundanceData, 'functionalOutput', functionalOutput, 'TCM', TCM, 'importanceValues', importanceValues, 'crossValidatedCoefficients', crossValidatedCoefficients);

% If usimg mock data, save synthetic coefficients and accuracy in addition
if ~usingRealData(settings)
    [results.syntheticCoefficients, results.accuracy] = deal(syntheticCoefficients, accuracy);
end
end