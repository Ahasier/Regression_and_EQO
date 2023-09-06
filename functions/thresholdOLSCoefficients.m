function [thresholedCoefficients, outSampleError] = thresholdOLSCoefficients(testData, testOutput, coefficients, threshold)
% initialize variables
weights = setWeightsRange();

% Compute out-of-sample mean squared error and r^2 on the test data to pick an
% optimal one with minimal MSE
[thresholedCoefficients, outSampleError] = pickOptimalWeightForThresholding(weights, threshold, coefficients, testData, testOutput);
end