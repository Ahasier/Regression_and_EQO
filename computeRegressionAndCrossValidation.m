function [allThresholdedCoefficients, avgBestCoefficients, allOutSampleErrors, avgBestOutSampleError, allOutSampleR2, allOptimalThresholds, coefficientsStdDev] = computeRegressionAndCrossValidation(abundanceData, functionalOutput, numPermutations, regressionMethod, Beta0)
for i = 1:numPermutations
    % Step 1: Split data
    [trainingData, testData, trainingOutput, testOutput] = splitData(abundanceData, functionalOutput);
    
    % Step 2: Regression on training data
    coefficients = computeRegression(trainingData, trainingOutput, regressionMethod, settings);
    inSampleError = computeSquaredError(trainingData, trainingOutput, coefficients); % For diagnostics
    
    % Step 3: Cross-validation on testing data
    [crossValidateThresholds, thresholdedCoefficients, bestCoefficients, outSampleErrors, bestOutSampleError, R2OutSamples, optimalThreshold] = computeCrossValidation(testData, testOutput, coefficients, regressionMethod, Beta0);
    
    % Step 4: Store results
    if i == 1 % If output variables are not exist, initialize them
        [allThresholdedCoefficients, allBestCoefficients, allOutSampleErrors, allBestOutSampleErrors, allOutSampleR2, allOptimalThresholds] = initializeCrossValidationResults(abundanceData, regressionMethod, coefficients, numPermutations, crossValidateThresholds);
    end
    [allThresholdedCoefficients, allBestCoefficients, allOutSampleErrors, allBestOutSampleErrors, allOutSampleR2, allOptimalThresholds] = storeCrossValidationResultsOverPermutations(i, regressionMethod, thresholdedCoefficients, bestCoefficients, outSampleErrors, bestOutSampleError, R2OutSamples, optimalThreshold, allThresholdedCoefficients, allBestCoefficients, allOutSampleErrors, allBestOutSampleErrors, allOutSampleR2, allOptimalThresholds);
end

% Aggregate results
avgBestCoefficients = mean(allBestCoefficients, 2);
coefficientsStdDev = std(allBestCoefficients, 0, 2);
avgBestOutSampleError = mean(allBestOutSampleErrors);
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

function [allThresholdedCoefficients, allBestCoefficients, allOutSampleErrors, allBestOutSampleErrors, allOutSampleR2, allOptimalThresholds] = storeCrossValidationResultsOverPermutations(idx, regressionMethod, thresholdedCoefficients, bestCoefficients, outSampleErrors, bestOutSampleError, R2OutSamples, optimalThreshold, allThresholdedCoefficients, allBestCoefficients, allOutSampleErrors, allBestOutSampleErrors, allOutSampleR2, allOptimalThresholds)
if strcmp(regressionMethod, 'OLS')
    allThresholdedCoefficients(:, :, idx) = thresholdedCoefficients;
    allBestCoefficients(:, idx) = bestCoefficients;
    allOutSampleErrors(:,idx) = outSampleErrors;
    allBestOutSampleErrors(idx) = bestOutSampleError;
    allOutSampleR2(:,idx) = R2OutSamples;
    allOptimalThresholds(idx) = optimalThreshold;
else
    allThresholdedCoefficients(:, :, :, idx) = thresholdedCoefficients;
    allBestCoefficients(:, :, idx) = bestCoefficients;
    allOutSampleErrors(:, :, idx) = outSampleErrors;
    allBestOutSampleErrors(:, idx) = bestOutSampleError;
    allOutSampleR2(:, :, idx) = R2OutSamples;
    allOptimalThresholds(:, idx) = optimalThreshold;
end
end

function [allThresholdedCoefficients, allBestCoefficients, allOutSampleErrors, allBestOutSampleErrors, allOutSampleR2, allOptimalThresholds] = initializeCrossValidationResults(abundanceData, regressionMethod, coefficients, numPermutations, crossValidateThresholds)
numThresholds = length(crossValidateThresholds);
if strcmp(regressionMethod, 'OLS')
    allThresholdedCoefficients = zeros(size(abundanceData,2), numThresholds, numPermutations);
    allBestCoefficients = zeros(size(abundanceData,2), numPermutations);
    allOutSampleErrors = zeros(numThresholds, numPermutations);
    allBestOutSampleErrors = zeros(1, numPermutations);
    allOutSampleR2 = zeros(numThresholds,numPermutations);
    allOptimalThresholds = zeros(1, numPermutations);
else
    numLambda = size(coefficients, 2);
    allThresholdedCoefficients = zeros(size(abundanceData,2), numLambda, numThresholds, numPermutations);
    allBestCoefficients = zeros(size(abundanceData,2), numLambda, numPermutations);
    allOutSampleErrors = zeros(numLambda, numThresholds, numPermutations);
    allBestOutSampleErrors = zeros(numLambda,  numPermutations);
    allOutSampleR2 = zeros(numLambda, numThresholds,numPermutations);
    allOptimalThresholds = zeros(numLambda,  numPermutations);
end
end