function [allThresholdedCoefficients, bestCoefficients, allOutSampleErrors, bestOutSampleError, allR2OutSample, optimalThreshold] = computeCrossValidation(crossValidateThresholds, testData, testOutput, estimatedCoefficients, regressionMethod)
% COMPUTECROSSVALIDATION performs cross-validation for the provided test data, estimated regression coefficients, and specified regression method.
% 
% INPUTS:
%   testData: Matrix containing the test data.
%   testOutput: Vector containing the corresponding test output.
%   estimatedCoefficients: Vector of regression coefficients estimated from the training phase.
%   regressionMethod: String specifying the regression method ('OLS', etc.).
%   beta0: Initial value for setting cross-validation thresholds.
% 
% OUTPUTS:
%   Various metrics related to cross-validation, including thresholded coefficients, errors, R-squared values, and the optimal threshold.

% Handle OLS regression method
if strcmp(regressionMethod, 'OLS')
    [allThresholdedCoefficients, allOutSampleErrors, allR2OutSample] = handleCrossValidationForOLS(crossValidateThresholds, testData, testOutput, estimatedCoefficients);
else
    % Handle other regression methods
    [allThresholdedCoefficients, allOutSampleErrors, allR2OutSample] = handleCrossValidationForOthers(crossValidateThresholds, testData, testOutput, estimatedCoefficients);
end
% Find the threshold that gives the minimum out-of-sample error
[bestOutSampleError, idx] = min(allOutSampleErrors(:));
bestCoefficients = allThresholdedCoefficients(:, idx);
optimalThreshold = findOptimalThreshold(regressionMethod, crossValidateThresholds, idx, estimatedCoefficients);
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

% Loop through all regresssion (lambda) parameters
for l = 1:len1
    % Get estimated coefficients at the current regresssion (lambda) parameters
    currentCoefficients = estimatedCoefficients(:, l);
    
    % For each threshold, compute out-of-sample error and r^2
    for k = 1:len2
        % Get the threshold at each iteration
        threshold = crossValidateThresholds(k);
        
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

function optimalThreshold = findOptimalThreshold(regressionMethod, crossValidateThresholds, idx, estimatedCoefficients)
if strcmp(regressionMethod, 'OLS')
    optimalThreshold = crossValidateThresholds(idx);
else
    lenLambda = size(estimatedCoefficients, 2);
    optimalThreshold = crossValidateThresholds(ceil(idx/lenLambda));
end
end