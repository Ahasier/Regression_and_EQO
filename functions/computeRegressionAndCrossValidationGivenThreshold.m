function [avgBestCoefficients, resultsForDiagnostics] = computeRegressionAndCrossValidationGivenThreshold(abundanceData, functionalOutput, numPermutations, regressionMethod, settings)
% Initialize output results variables
allCoefficients = zeros(size(abundanceData,2), numPermutations);
allOutSampleErrors = zeros(1, numPermutations);

% Loop over various permutations
for i = 1:numPermutations
    % Step 1: Split data
    [trainingData, testData, trainingOutput, testOutput] = splitData(abundanceData, functionalOutput);
    
    % Step 2: Regression on training data
    coefficients = computeRegression(trainingData, trainingOutput, regressionMethod, settings);
    inSampleError = computeSquaredError(trainingData, trainingOutput, coefficients); % For diagnostics
    
    % Step 3: Cross-validation on testing data
    % Handle OLS regression method
    if strcmp(regressionMethod, 'OLS')
        [thresholedCoefficients, outSampleError] = thresholdOLSCoefficients(testData, testOutput, coefficients, settings.Threshold);
    else
        % Handle other regression methods
        [thresholedCoefficients, ~, outSampleError] = crossValidationForLASSO(testData, testOutput, coefficients, settings.maxLambda, settings.Threshold);
        
    end
    
    % Step 4: Store results
    [allCoefficients(:, i), allOutSampleErrors(i)] = storeThresholdedResults(thresholedCoefficients, outSampleError);
end

% Aggregate results
avgBestCoefficients = median(allCoefficients, 2);
avgBestOutSampleError = mean(allOutSampleErrors);

% Store results as a structure
resultsForDiagnostics = struct('allCoefficients', allCoefficients, 'allErrorsAtAllThresholds', allOutSampleErrors, 'MeanSquaredErrorOutOfSample', avgBestOutSampleError);
end

% Helper functions
function [thisCoefficient, thisOutSampleErrors] = storeThresholdedResults(thresholedCoefficients, outSampleError)
thisCoefficient = thresholedCoefficients;
thisOutSampleErrors = outSampleError;
end

function [thresholedCoefficients, outSampleError] = thresholdOLSCoefficients(testData, testOutput, coefficients, threshold)
thresholedCoefficients = coefficients > threshold;
outSampleError = computeSquaredError(testData, testOutput, thresholedCoefficients);
end