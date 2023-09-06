function [optimallyWeightedThresholdedCoefficients, optimallyWeightedThresholdedOutSampleError] = pickOptimalWeightForThresholding(weights, threshold, estimatedCoefficients, testData, testOutput)
% Initialize weighted coefficients and out of sample mean squared errors over different weights
lenWeights = length(weights);
numTaxa = size(estimatedCoefficients, 1);
OutSampleErrorsOverDifferentWeights = zeros(1,lenWeights);
thresholdedCoefficientsOverDifferentWeights = zeros(numTaxa, lenWeights);

% Compute the binary vector of whether coefficients exceed threshold
binaryCoefficients = estimatedCoefficients >= threshold;
thresholdedCoefficients = binaryCoefficients.*estimatedCoefficients;

% Loop over different weights, threshold coefficients and compute out-of-sample
% errors on the test data
for w = 1:lenWeights
    weight = weights(w);
    weightedThresholdedCoefficients = weight.*thresholdedCoefficients;
    OutSampleErrorsOverDifferentWeights(w) = computeSquaredError(testData, testOutput, weightedThresholdedCoefficients);
    thresholdedCoefficientsOverDifferentWeights(:, w) = weightedThresholdedCoefficients;
end

% Pick the optimal weight
[optimallyWeightedThresholdedOutSampleError, idxOptimalWeight] = min(OutSampleErrorsOverDifferentWeights);

% The out-of-sample error and r^2 are then those cumputed from the optimally
% weighted coefficients.
optimallyWeightedThresholdedCoefficients = thresholdedCoefficientsOverDifferentWeights(:, idxOptimalWeight);
end