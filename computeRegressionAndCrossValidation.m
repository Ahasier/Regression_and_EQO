function [avgBestCoefficients, resultsForDiagnostics] = computeRegressionAndCrossValidation(abundanceData, functionalOutput, numPermutations, regressionMethod, settings)
% COMPUTEREGRESSIONANDCROSSVALIDATION integrates the regression and cross-validation processes.
% 
% INPUTS:
%   abundanceData: Matrix containing the abundance data of different taxa across samples.
%   functionalOutput: Vector containing the specific functional output in the samples.
%   numPermutations: Number of times to perform the integrated process.
%   regressionMethod: String specifying the regression method ('L0', 'OLS', etc.).
%   settings: A structure with additional parameters and settings.
%
% OUTPUTS:
%   Various metrics related to the integrated regression and cross-validation process, including coefficients, out-of-sample errors, R-squared values, and optimal thresholds.

% Set threshold values for cross validation to go over
crossValidateThresholds = setCrossValidationThreshold(beta0);

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
    [thresholdedCoefficients, bestThresholedCoefficients, outSampleErrors, bestOutSampleError, R2OutSamples, optimalThreshold] = computeCrossValidation(crossValidateThresholds, testData, testOutput, coefficients, regressionMethod, settings.Beta0);
    
    % Step 4: Store results
    if strcmp(regressionMethod,'OLS')
        [allCoefficients(:, :, idx), allBestCoefficients(:, idx), allOutSampleErrors(:,idx), allBestOutSampleErrors(idx), allOutSampleR2(:,idx), allOptimalThresholds(idx)] = storeCrossValidationResults(thresholdedCoefficients, bestThresholedCoefficients, outSampleErrors, bestOutSampleError, R2OutSamples, optimalThreshold);
    else
        [allCoefficients(:, :, :, idx), allBestCoefficients(:, idx), allOutSampleErrors(:, :, idx), allBestOutSampleErrors(idx), allOutSampleR2(:, :, idx), allOptimalThresholds(:, idx)] = storeCrossValidationResults(thresholdedCoefficients, bestThresholedCoefficients, outSampleErrors, bestOutSampleError, R2OutSamples, optimalThreshold);
    end
end

% Aggregate results
avgBestCoefficients = median(allBestCoefficients, 2);
coefficientsStdDev = std(allBestCoefficients, 0, 2);
avgBestOutSampleError = mean(allBestOutSampleErrors);

% Store results as a structure
resultsForDiagnostics = struct('allCoefficients', allCoefficients, 'allBestCoefficients', allBestCoefficients, 'allErrorsAtAllThresholds', allOutSampleErrors, 'allOutSampleR2AtAllThresholds', allOutSampleR2, 'allOptimalThresholds', allOptimalThresholds, 'coefficientsStdDev', coefficientsStdDev, 'MeanSquaredErrorOutOfSample', avgBestOutSampleError);
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

function crossValidateThresholds = setCrossValidationThreshold(beta0)
crossValidateThresholds = 0.1:0.1:2*beta0;
end