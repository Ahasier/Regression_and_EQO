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