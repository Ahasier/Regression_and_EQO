% Perform cross-validation for regressions
function [allThresholdedCoefficients, allOutSampleErrors, allR2OutSample] = crossValidateOverThresholds(crossValidateThresholds, testData, testOutput, estimatedCoefficients)
% initialize variables
len = length(crossValidateThresholds);
allOutSampleErrors = zeros(1,len);
allR2OutSample = zeros(1,len);
allThresholdedCoefficients = logical(zeros(size(estimatedCoefficients, 1), len));
% weights = setWeightsRange();

% For each threshold, compute out-of-sample error
for l = 1:len
    % Get the threshold at each iteration
    threshold = crossValidateThresholds(l);
    
    % Threshold coefficients and multiple by different weights, compute out-of-sample
    % error on the test data to pick an optimal one.
    [binaryCoefficients, weightedThresholdedCoefficients, outSampleError] = weightForThresholding(threshold, estimatedCoefficients, testData, testOutput);
    
    % Store results
    allR2OutSample(l) = computeRSquared(testData, testOutput, weightedThresholdedCoefficients);
    allOutSampleErrors(l) = outSampleError;
    allThresholdedCoefficients(:, l) = binaryCoefficients;
end
end