function [crossValidateThresholds, allCoefficients, bestCoefficients, allOutSampleErrors, bestOutSampleError, allR2OutSample, optimalThreshold] = computeCrossValidation(testData, testOutput, estimatedCoefficients, regressionMethod, beta0)
% Handle OLS regression method
if strcmp(regressionMethod, 'OLS')
    [crossValidateThresholds, allCoefficients, allOutSampleErrors, allR2OutSample] = handleCrossValidationForOLS(testData, testOutput, estimatedCoefficients, beta0);
else
    % Handle other regression methods
    [crossValidateThresholds, allCoefficients, allOutSampleErrors, allR2OutSample] = handleCrossValidationForOthers(testData, testOutput, estimatedCoefficients, beta0);
end

% Find the threshold that gives the minimum out-of-sample error
[bestOutSampleError, idx] = min(allOutSampleErrors);
bestCoefficients = allCoefficients(:, idx);
optimalThreshold = crossValidateThresholds(idx);
end

%% Helper functions
% Perform cross-validation for Ordinary Least Squares regression
function [crossValidateThresholds, allThresholdedCoefficients, allOutSampleErrors, allR2OutSample] = handleCrossValidationForOLS(testData, testOutput, estimatedCoefficients, beta0)
% initialize variables
crossValidateThresholds = 0:0.1:2*beta0;
numTaxa = size(estimatedCoefficients, 1);
len = length(crossValidateThresholds);
allOutSampleErrors = zeros(1,len);
allR2OutSample = zeros(1,len);
allThresholdedCoefficients = zeros(size(estimatedCoefficients, 1), len);
weights = 0:0.01:2;
lenWeights = length(weights);

% For each threshold, compute out-of-sample error
for l = 1:len
    % Initialize weighted coefficients and out of sample mean squared errors over different weights
    OutSampleErrorsOverDifferentWeights = zeros(1,lenWeights);
    thresholdedCoefficientsOverDifferentWeights = zeros(numTaxa, lenWeights);
    
    threshold = crossValidateThresholds(l);
    binaryCoefficients = estimatedCoefficients >= threshold;
    % Threshold coefficients and multiple by different weights, compute out-of-sample 
    % error on the test data to pick an optimal one.
    for w = 1:lenWeights
        weight = weights(w);
        thresholdedCoefficients = weight.*binaryCoefficients.*estimatedCoefficients;
        OutSampleErrorsOverDifferentWeights(w) = computeSquaredError(testData, testOutput, thresholdedCoefficients);
        thresholdedCoefficientsOverDifferentWeights(:, w) = thresholdedCoefficients;
    end
    
    % Pick the optimal weight
    [~, idxOptimalWeight] = min(OutSampleErrorsOverDifferentWeights);
    
    % The out-of-sample error and r^2 are then those cumputed from the optimally
    % weighted coefficients.
    optimallyWeightedThresholdedCoefficients = thresholdedCoefficientsOverDifferentWeights(:, idxOptimalWeight);
    allR2OutSample(l) = computeRSquared(testData, testOutput, optimallyWeightedThresholdedCoefficients);
    allOutSampleErrors(l) = OutSampleErrorsOverDifferentWeights(idxOptimalWeight);
    allThresholdedCoefficients(:, l) = optimallyWeightedThresholdedCoefficients;
end
end

% Handle other regression methods for cross-validation
function [crossValidateThresholds, allThresholdedCoefficients, allOutSampleErrors, allR2OutSample] = handleCrossValidationForOthers(testData, testOutput, estimatedCoefficients, beta0)
% initialize variables
crossValidateThresholds = 0:0.1:2*beta0;
numTaxa = size(estimatedCoefficients, 1);
len1 = size(estimatedCoefficients, 2);
len2 = length(crossValidateThresholds);
allOutSampleErrors = zeros(len1, len2);
allR2OutSample = zeros(len1, len2);
allThresholdedCoefficients = zeros(size(estimatedCoefficients, 1), len1, len2);
weights = 0.9:0.001:1.1;
lenWeights = length(weights);

% Iterate through all regresssion (lambda) parameters
for l = 1:len1
    % For each threshold, compute out-of-sample error and r^2
    for k = 1:len2
        % Initialize weighted coefficients and out of sample mean squared errors over different weights
        OutSampleErrorsOverDifferentWeights = zeros(1,lenWeights);
        thresholdedCoefficientsOverDifferentWeights = zeros(numTaxa, lenWeights);
        
        threshold = crossValidateThresholds(l);
        binaryCoefficients = estimatedCoefficients(:, k) >= threshold;
        % Threshold coefficients and multiple by different weights, compute out-of-sample
        % error on the test data to pick an optimal one.
        for w = 1:lenWeights
            weight = weights(w);
            thresholdedCoefficients = weight.*binaryCoefficients.*estimatedCoefficients;
            OutSampleErrorsOverDifferentWeights(w) = computeSquaredError(testData, testOutput, thresholdedCoefficients);
            thresholdedCoefficientsOverDifferentWeights(:, w) = thresholdedCoefficients;
        end
        % Pick the optimal weight
        [~, idxOptimalWeight] = min(OutSampleErrorsOverDifferentWeights);
        
        % The out-of-sample error and r^2 are then those cumputed from the optimally
        % weighted coefficients.
        optimallyWeightedThresholdedCoefficients = thresholdedCoefficientsOverDifferentWeights(:, idxOptimalWeight);
        allR2OutSample(l, k) = computeRSquared(testData, testOutput, optimallyWeightedThresholdedCoefficients);
        allOutSampleErrors(l, k) = OutSampleErrorsOverDifferentWeights(idxOptimalWeight);
        allThresholdedCoefficients(:, l, k) = optimallyWeightedThresholdedCoefficients;
    end
end
end