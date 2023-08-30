function [avgCoefficients, coefficientsStdDev, avgOutSampleError] = computeRegressionAndCrossValidation(abundanceData, functionalOutput, numPermutations, regressionMethod, Beta0)
allCoefficients = [];
allErrors = [];
for i = 1:numPermutations
    % Step 1: Split data
    [trainingData, testData, trainingOutput, testOutput] = splitData(abundanceData, functionalOutput);
    
    % Step 2: Regression on training data
    coefficients = computeRegression(trainingData, trainingOutput, regressionMethod, Beta0);
    inSampleError = computeError(trainingData, trainingOutput, coefficients); % For diagnostics
    
    % Step 3: Cross-validation on testing data
    [thresholdedCoefficients, outSampleError] = computeCrossValidation(testData, testOutput, coefficients);
    
    % Step 4: Store results
    allCoefficients = [allCoefficients; thresholdedCoefficients];
    allErrors = [allErrors; outSampleError];
end

% Aggregate results
avgCoefficients = mean(allCoefficients, 1);
coefficientsStdDev = std(allCoefficients, 0, 1);
avgOutSampleError = mean(allErrors);
end