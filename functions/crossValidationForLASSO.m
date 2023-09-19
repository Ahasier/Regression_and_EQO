function [optimalCoefficients, optimalLambda, bestOutSampleError] = crossValidationForLASSO(testData, testOutput, coefficients, maxLambda)
% initialize variables
Lambda = setLambdaRange(maxLambda);
len1 = length(Lambda);
allOutSampleErrors = zeros(len1, 1);
allR2OutSample = zeros(len1, 1);
allCoefficients = zeros(size(coefficients, 1), len1);

% Loop through all regresssion (lambda) parameters
for l = 1:len1
    % Get estimated coefficients at the current regresssion (lambda) parameters
    currentCoefficients = coefficients(:, l);
        
    % Store results
    allR2OutSample(l) = computeRSquared(testData, testOutput, currentCoefficients);
    allOutSampleErrors(l) = computeSquaredError(testData, testOutput, currentCoefficients);
    allCoefficients(:, l) = currentCoefficients;
end
[bestOutSampleError, idx] = min(allOutSampleErrors(:));
optimalCoefficients = allCoefficients(:, idx);
optimalLambda = Lambda(idx);
end