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

if isfield(settings, 'RealAbd') && strcmp(settings.RealAbd, 'On')
    results = processData(settings, numPermutations, noiseLevel, numSamples, regressionMethod);
else
    % Load tree data
    treeData = loadTreeData();
    
    % Process data to get results
    results = processData(settings, numPermutations, noiseLevel, numSamples, regressionMethod, treeData, phylogenyDependency, numberOfTaxaInAGroup);
end

% Save results to disk
saveResults(results, regressionMethod, fullIdentifier, numberOfTaxaInAGroup, numSamples, settings, meshGrid);
end

%% Helper functions
% Load tree data
function treeData = loadTreeData()
global paths

load([paths.data,'tree100taxa.mat'], 'tree100');
treeData = tree100;
end

function results = processData(settings, numPermutations, noiseLevel, numSamples, regressionMethod, varargin)
% Check the number of inputs
numInputs = length(varargin);

% 1. Generate Data
if numInputs == 3
    treeData = varargin{1};
    phylogenyDependency = varargin{2};
    numberOfTaxaInAGroup = varargin{3};
    [abundanceData, functionalOutput, varargout] = generateData(treeData, phylogenyDependency, numberOfTaxaInAGroup, noiseLevel, numSamples, settings);
elseif numInputs == 0
    index = settings.index;
    [abundanceData, functionalOutput, varargout] = generateDataFromRealSample(index, noiseLevel, settings, numSamples);
else
    error("invalid number of varargin.")
end

% If usimg mock data, retrieve the synthetic coefficients from varargout
if ~usingRealData(settings)
    syntheticCoefficients = varargout{1};
    if numInputs == 0
        numberOfTaxaInAGroup = sum(syntheticCoefficients > 0);
    end
end

% If using incorporating phylogenetic feature, retrieve the grouped
% abundance and grouping indices, and pass them to
% `computeRegressionAndCrossValidation`
% If using EQO, use `numberOfTaxaInAGroup` to set an appropriate value for
% N_max, to reduce the computation time. (The two conditions don't coexist)
if useExtraFeatures(settings)
    extraPhyloVars = varargout{2};
    varargin = extraPhyloVars;
elseif strcmp(regressionMethod, 'EQO')
    varargin = numberOfTaxaInAGroup;
else
    varargin = {};
end

% 2. Compute regression and cross-validation over different thresholds on data
[TCM, importanceValues, crossValidatedCoefficients] = computeRegressionAndCrossValidation(abundanceData, functionalOutput, numPermutations, regressionMethod, settings, varargin);

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