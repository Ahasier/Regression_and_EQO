function [optimalCoefficients, optimalLambda, bestOutSampleError] = crossValidationForLASSO(testData, testOutput, coefficients, maxLambda, threshold)
% initialize variables
Lambda = setLambdaRange(maxLambda);
len1 = length(Lambda);
allOutSampleErrors = zeros(len1, 1);
allR2OutSample = zeros(len1, 1);
allCoefficients = zeros(size(coefficients, 1), len1);
weights = setWeightsRange();

% Loop through all regresssion (lambda) parameters
for l = 1:len1
    % Get estimated coefficients at the current regresssion (lambda) parameters
    currentCoefficients = coefficients(:, l);
    
    % Compute out-of-sample mean squared error and r^2 on the test data to pick an
    % optimal one with minimal MSE
    [optimallyWeightedCoefficients, optimallyWeightedOutSampleError] = pickOptimalWeightForThresholding(weights, threshold, currentCoefficients, testData, testOutput);
    
    % Store results
    allR2OutSample(l) = computeRSquared(testData, testOutput, optimallyWeightedCoefficients);
    allOutSampleErrors(l) = optimallyWeightedOutSampleError;
    allCoefficients(:, l) = optimallyWeightedCoefficients;
end
[bestOutSampleError, idx] = min(allOutSampleErrors(:));
optimalCoefficients = allCoefficients(:, idx);
optimalLambda = Lambda(idx);
end