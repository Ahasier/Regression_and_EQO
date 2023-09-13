function [avgBestCoefficients, resultsForDiagnostics] = computeRegressionAndCrossValidationViaThresholding(abundanceData, functionalOutput, numPermutations, regressionMethod, settings)
% Set threshold values for cross validation to go over
crossValidateThresholds = setCrossValidationThreshold(settings.Beta0);

% Initialize output results variables
[allCoefficients, allBestCoefficients, allOutSampleErrors, allBestOutSampleErrors, allOutSampleR2, allOptimalThresholds] = initializeCrossValidationResults(abundanceData, regressionMethod, numPermutations, crossValidateThresholds, settings);

% Loop over various permutations
for i = 1:numPermutations
    % Step 1: Split data
    [trainingData, testData, trainingOutput, testOutput] = splitData(abundanceData, functionalOutput);
    
    % Step 2: Regression on training data
    coefficients = computeRegression(trainingData, trainingOutput, regressionMethod, settings);
    inSampleError = computeSquaredError(trainingData, trainingOutput, coefficients); % For diagnostics
    
    % Step 3: Cross-validation on testing data
    [thresholdedCoefficients, bestThresholedCoefficients, outSampleErrors, bestOutSampleError, R2OutSamples, optimalThreshold] = computeCrossValidation(crossValidateThresholds, testData, testOutput, coefficients, regressionMethod);
    
    % Step 4: Store results
    if strcmp(regressionMethod,'OLS')
        [allCoefficients(:, :, i), allBestCoefficients(:, i), allOutSampleErrors(:, i), allBestOutSampleErrors(i), allOutSampleR2(:, i), allOptimalThresholds(i)] = storeCrossValidationResults(thresholdedCoefficients, bestThresholedCoefficients, outSampleErrors, bestOutSampleError, R2OutSamples, optimalThreshold);
    else
        [allCoefficients(:, :, :, i), allBestCoefficients(:, i), allOutSampleErrors(:, :, i), allBestOutSampleErrors(i), allOutSampleR2(:, :, i), allOptimalThresholds(:, i)] = storeCrossValidationResults(thresholdedCoefficients, bestThresholedCoefficients, outSampleErrors, bestOutSampleError, R2OutSamples, optimalThreshold);
    end
end

% Aggregate results and store them as a structure
[avgBestCoefficients, resultsForDiagnostics] = aggregateCrossValidationResults(allCoefficients, allOutSampleErrors, allOutSampleR2, allOptimalThresholds, allBestCoefficients, allBestOutSampleErrors);
end

%% Helper functions
function [thisCoefficients, thisBestCoefficients, thisOutSampleErrors, thisBestOutSampleErrors, thisOutSampleR2, thisOptimalThresholds] = storeCrossValidationResults(coefficients, bestCoefficients, outSampleErrors, bestOutSampleError, R2OutSamples, optimalThreshold)
thisCoefficients = coefficients;
thisBestCoefficients = bestCoefficients;
thisOutSampleErrors = outSampleErrors;
thisBestOutSampleErrors = bestOutSampleError;
thisOutSampleR2 = R2OutSamples;
thisOptimalThresholds = optimalThreshold;
end

function [allCoefficients, allBestCoefficients, allOutSampleErrors, allBestOutSampleErrors, allOutSampleR2, allOptimalThresholds] = initializeCrossValidationResults(abundanceData, regressionMethod, numPermutations, crossValidateThresholds, settings)
numThresholds = length(crossValidateThresholds);
if strcmp(regressionMethod, 'OLS')
    allCoefficients = zeros(size(abundanceData,2), numThresholds, numPermutations);
    allBestCoefficients = zeros(size(abundanceData,2), numPermutations);
    allOutSampleErrors = zeros(numThresholds, numPermutations);
    allBestOutSampleErrors = zeros(1, numPermutations);
    allOutSampleR2 = zeros(numThresholds,numPermutations);
    allOptimalThresholds = zeros(1, numPermutations);
else
    Lambda = setLambdaRange(settings.maxLambda);
    numLambda = length(Lambda);
    allCoefficients = zeros(size(abundanceData,2), numLambda, numThresholds, numPermutations);
    allBestCoefficients = zeros(size(abundanceData,2), numPermutations);
    allOutSampleErrors = zeros(numLambda, numThresholds, numPermutations);
    allBestOutSampleErrors = zeros(1,  numPermutations);
    allOutSampleR2 = zeros(numLambda, numThresholds,numPermutations);
    allOptimalThresholds = zeros(numLambda,  numPermutations);
end
end

function [avgBestCoefficients, resultsForDiagnostics] = aggregateCrossValidationResults(allCoefficients, allOutSampleErrors, allOutSampleR2, allOptimalThresholds, allBestCoefficients, allBestOutSampleErrors)
% Aggregate results
avgBestCoefficients = median(allBestCoefficients, 2);
coefficientsStdDev = std(allBestCoefficients, 0, 2);
avgBestOutSampleError = mean(allBestOutSampleErrors);

% Store results as a structure
resultsForDiagnostics = struct('allCoefficients', allCoefficients, 'allBestCoefficients', allBestCoefficients, 'allErrorsAtAllThresholds', allOutSampleErrors, 'allOutSampleR2AtAllThresholds', allOutSampleR2, 'allOptimalThresholds', allOptimalThresholds, 'coefficientsStdDev', coefficientsStdDev, 'MeanSquaredErrorOutOfSample', avgBestOutSampleError);
end