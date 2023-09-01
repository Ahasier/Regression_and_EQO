function [crossValidateThresholds, allThresholdedCoefficients, bestCoefficients, allOutSampleErrors, bestOutSampleError, allR2OutSample, optimalThreshold] = computeCrossValidation(testData, testOutput, estimatedCoefficients, regressionMethod, beta0)
% Set threshold values for cross validation to go over
crossValidateThresholds = setCrossValidationThreshold(beta0);

% Handle OLS regression method
if strcmp(regressionMethod, 'OLS')
    [allThresholdedCoefficients, allOutSampleErrors, allR2OutSample] = handleCrossValidationForOLS(crossValidateThresholds, testData, testOutput, estimatedCoefficients);
else
    % Handle other regression methods
    [allThresholdedCoefficients, allOutSampleErrors, allR2OutSample] = handleCrossValidationForOthers(crossValidateThresholds, testData, testOutput, estimatedCoefficients);
end

% Find the threshold that gives the minimum out-of-sample error
[bestOutSampleError, idx] = min(allOutSampleErrors);
bestCoefficients = allThresholdedCoefficients(:, idx);
optimalThreshold = crossValidateThresholds(idx);
end

%% Helper functions
% Perform cross-validation for Ordinary Least Squares regression
function [allThresholdedCoefficients, allOutSampleErrors, allR2OutSample] = handleCrossValidationForOLS(crossValidateThresholds, testData, testOutput, estimatedCoefficients)
% initialize variables
len = length(crossValidateThresholds);
allOutSampleErrors = zeros(1,len);
allR2OutSample = zeros(1,len);
allThresholdedCoefficients = zeros(size(estimatedCoefficients, 1), len);
weights = setWeightsRange();

% For each threshold, compute out-of-sample error
for l = 1:len
    % Get the threshold at each iteration
    threshold = crossValidateThresholds(l);
    
    % Threshold coefficients and multiple by different weights, compute out-of-sample
    % error on the test data to pick an optimal one.
    [optimallyWeightedThresholdedCoefficients, optimallyWeightedThresholdedOutSampleError] = pickOptimalWeightForThresholding(weights, threshold, estimatedCoefficients, testData, testOutput);
    
    % Store results
    allR2OutSample(l) = computeRSquared(testData, testOutput, optimallyWeightedThresholdedCoefficients);
    allOutSampleErrors(l) = optimallyWeightedThresholdedOutSampleError;
    allThresholdedCoefficients(:, l) = optimallyWeightedThresholdedCoefficients;
end
end

% Handle other regression methods for cross-validation
function [allThresholdedCoefficients, allOutSampleErrors, allR2OutSample] = handleCrossValidationForOthers(crossValidateThresholds, testData, testOutput, estimatedCoefficients)
% initialize variables
len1 = size(estimatedCoefficients, 2);
len2 = length(crossValidateThresholds);
allOutSampleErrors = zeros(len1, len2);
allR2OutSample = zeros(len1, len2);
allThresholdedCoefficients = zeros(size(estimatedCoefficients, 1), len1, len2);
weights = setWeightsRange();

% Iterate through all regresssion (lambda) parameters
for l = 1:len1
    % Get estimated coefficients at the current regresssion (lambda) parameters
    currentCoefficients = estimatedCoefficients(:, l);
    
    % For each threshold, compute out-of-sample error and r^2
    for k = 1:len2
        % Get the threshold at each iteration
        threshold = crossValidateThresholds(l);
        
        % Threshold coefficients and multiple by different weights, compute out-of-sample
        % error on the test data to pick an optimal one.
        [optimallyWeightedThresholdedCoefficients, optimallyWeightedThresholdedOutSampleError] = pickOptimalWeightForThresholding(weights, threshold, currentCoefficients, testData, testOutput);
        
        % Store results
        allR2OutSample(l, k) = computeRSquared(testData, testOutput, optimallyWeightedThresholdedCoefficients);
        allOutSampleErrors(l, k) = optimallyWeightedThresholdedOutSampleError;
        allThresholdedCoefficients(:, l, k) = optimallyWeightedThresholdedCoefficients;
    end
end
end

function [optimallyWeightedThresholdedCoefficients, optimallyWeightedThresholdedOutSampleError] = pickOptimalWeightForThresholding(weights, threshold, estimatedCoefficients, testData, testOutput)
% Initialize weighted coefficients and out of sample mean squared errors over different weights
lenWeights = length(weights);
numTaxa = size(estimatedCoefficients, 1);
OutSampleErrorsOverDifferentWeights = zeros(1,lenWeights);
thresholdedCoefficientsOverDifferentWeights = zeros(numTaxa, lenWeights);

% Compute the binary vector of whether coefficients exceed threshold
binaryCoefficients = estimatedCoefficients >= threshold;

% Iterate over different weights, threshold coefficients and compute out-of-sample
% errors on the test data
for w = 1:lenWeights
    weight = weights(w);
    thresholdedCoefficients = weight.*binaryCoefficients.*estimatedCoefficients;
    OutSampleErrorsOverDifferentWeights(w) = computeSquaredError(testData, testOutput, thresholdedCoefficients);
    thresholdedCoefficientsOverDifferentWeights(:, w) = thresholdedCoefficients;
end

% Pick the optimal weight
[optimallyWeightedThresholdedOutSampleError, idxOptimalWeight] = min(OutSampleErrorsOverDifferentWeights);

% The out-of-sample error and r^2 are then those cumputed from the optimally
% weighted coefficients.
optimallyWeightedThresholdedCoefficients = thresholdedCoefficientsOverDifferentWeights(:, idxOptimalWeight);
end

function crossValidateThresholds = setCrossValidationThreshold(beta0)
crossValidateThresholds = 0.1:0.1:2*beta0;
end

function weights = setWeightsRange()
weights = 0:0.001:20;
end