function [binaryCoefficients, weightedThresholdedCoefficients, outSampleError] = weightForThresholding(threshold, estimatedCoefficients, testData, testOutput)
% Compute the binary vector of whether coefficients exceed threshold
binaryCoefficients = estimatedCoefficients >= threshold;
thresholdedCoefficients = binaryCoefficients.*estimatedCoefficients;

% Find the optimal weight through a regression
P = testData * thresholdedCoefficients;

% The weighted predicted outpus are `y_pred = b * P`, where the values 
% of b are inferred from a regression.
weight = P \ testOutput;

% Get the weighted thresholded coefficients and the out of sample error
weightedThresholdedCoefficients = weight * thresholdedCoefficients;
outSampleError = computeSquaredError(testData, testOutput, weightedThresholdedCoefficients);
end