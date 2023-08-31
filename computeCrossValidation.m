function [allThresholdedCoefficients, allOutSampleErrors, bestCoefficients, bestOutSampleError, optimalThreshold] = computeCrossValidation(testData, testOutput, estimatedCoefficients, regressionMethod, beta0)
% Handle OLS regression method
if strcmp(regressionMethod, 'OLS')
    [crossValidateThresholds, allThresholdedCoefficients, allOutSampleErrors] = handleCrossValidationForOLS(testData, testOutput, estimatedCoefficients, beta0);
else
    % Handle other regression methods
    [crossValidateThresholds, allThresholdedCoefficients, allOutSampleErrors] = handleCrossValidationForOthers(testData, testOutput, estimatedCoefficients, beta0);
end

% Find the threshold that gives the minimum out-of-sample error
[bestOutSampleError, idx] = min(allOutSampleErrors);
bestCoefficients = allThresholdedCoefficients(:, idx);
optimalThreshold = crossValidateThresholds(idx);
end

%% Helper functions
% Perform cross-validation for Ordinary Least Squares regression
function [crossValidateThresholds, allThresholdedCoefficients, allOutSampleErrors, allR2OutSample] = handleCrossValidationForOLS(testData, testOutput, estimatedCoefficients, beta0)
% initialize variables
crossValidateThresholds = 0:0.1:2*beta0;
len = length(crossValidateThresholds);
allOutSampleErrors = zeros(1,len);
allR2OutSample = zeros(1,len);
allThresholdedCoefficients = zeros(size(estimatedCoefficients, 1), len);

% For each threshold, compute out-of-sample error
for l = 1:len
    threshold = crossValidateThresholds(l);
    binaryCoefficients = estimatedCoefficients >= threshold;
    thresholdedCoefficients = binaryCoefficients.*estimatedCoefficients;
    
    % Compute out-of-sample error and r^2 on the test data
    allR2OutSample(l) = computeRSquared(testData, testOutput, thresholdedCoefficients);
    allOutSampleErrors(l) = computeRSquared(testData, testOutput, thresholdedCoefficients);
    allThresholdedCoefficients(:, l) = thresholdedCoefficients;
end
end

% Handle other regression methods for cross-validation
function [crossValidateThresholds, allThresholdedCoefficients, allOutSampleErrors, allR2OutSample] = handleCrossValidationForOthers(testData, testOutput, estimatedCoefficients, beta0)
% initialize variables
crossValidateThresholds = 0:0.1:2*beta0;
len1 = size(estimatedCoefficients, 2);
len2 = length(crossValidateThresholds);
allOutSampleErrors = zeros(len1, len2);
allR2OutSample = zeros(len1, len2);
allThresholdedCoefficients = zeros(size(estimatedCoefficients, 1), len1, len2);

% Iterate through all regresssion (lambda) parameters
for l = 1:len1
    % For each threshold, compute out-of-sample error and r^2
    for k = 1:len2
        threshold = crossValidateThresholds(l);
        binaryCoefficients = estimatedCoefficients(:, k) >= threshold;
        thresholdedCoefficients = binaryCoefficients.*estimatedCoefficients(:, k);
        
        % Compute out-of-sample error on the test data
        allR2OutSample(l, k) = computeRSquared(testData, testOutput, thresholdedCoefficients);
        allOutSampleErrors(l, k) = computeSquaredError(testData, testOutput, thresholdedCoefficients);
        allThresholdedCoefficients(:, l, k) = thresholdedCoefficients;
    end
end
end