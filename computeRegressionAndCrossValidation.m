function [allCoefficientsAtAllThresholds, avgCoefficients, optimalThreshold, coefficientsStdDev, avgOutSampleError] = computeRegressionAndCrossValidation(abundanceData, functionalOutput, numPermutations, regressionMethod, Beta0)
allCoefficientsAtAllThresholds = [];
allErrorsAtAllThresholds = [];
allCoefficients = [];
allErrors = [];
for i = 1:numPermutations
    % Step 1: Split data
    [trainingData, testData, trainingOutput, testOutput] = splitData(abundanceData, functionalOutput);
    
    % Step 2: Regression on training data
    coefficients = computeRegression(trainingData, trainingOutput, regressionMethod, settings);
    inSampleError = computeSquaredError(trainingData, trainingOutput, coefficients); % For diagnostics
    
    % Step 3: Cross-validation on testing data
    [coefficientsAtAllThresholds, outSampleErrorsAtAllThresholds, thresholdedCoefficients, outSampleError, optimalThreshold] = computeCrossValidation(testData, testOutput, coefficients, regressionMethod, Beta0);
    
    % Step 4: Store results
    allCoefficientsAtAllThresholds(:, :, end + 1) = coefficientsAtAllThresholds;
    allCoefficients = [allCoefficients, thresholdedCoefficients];
    allErrorsAtAllThresholds(:, end + 1) = outSampleErrorsAtAllThresholds;
    allErrors = [allErrors; outSampleError];
end

% Aggregate results
avgCoefficients = mean(allCoefficients, 2);
coefficientsStdDev = std(allCoefficients, 0, 2);
avgOutSampleError = mean(allErrors);
end

%% Helper functions
function [trainingData, testData, trainingOutput, testOutput] = splitData(abundanceData, functionalOutput)
totalSamples = size(abundanceData, 1);
idx = randperm(totalSamples);

% Assuming a 50-50 split for training and testing
splitPoint = floor(0.5 * totalSamples);

trainingIdx = idx(1:splitPoint);
testIdx = idx(splitPoint+1:end);

trainingData = abundanceData(trainingIdx, :);
testData = abundanceData(testIdx, :);

trainingOutput = functionalOutput(trainingIdx);
testOutput = functionalOutput(testIdx);
end