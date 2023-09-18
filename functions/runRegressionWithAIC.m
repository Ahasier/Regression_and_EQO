function [binaryCoefficients, optimalGroupSize] = runRegressionWithAIC(trainingData, trainingOutput, regressionMethod, settings, varargin)
% Step (b) 1: Regression on training data
coefficients = computeRegression(trainingData, trainingOutput, regressionMethod, settings);

% If using regressions other than OLS (i.e, LASSO), cross-validate to get the optimal lambda
if ~strcmp(regressionMethod, 'OLS')
    testData = varargin{1};
    testOutput = varargin{2};
    [coefficients, ~] = crossValidationForLASSO(testData, testOutput, coefficients, settings.maxLambda);
end

% Step (b) 2: Calculate the AIC values under different group sizes to find the optimal one
aicValues = evaluateAIC(trainingData, trainingOutput, coefficients);
[~, optimalGroupSize] = findMinimalAic(aicValues);

% Step (b) 3: Calculate the binary results sparsified from AIC step
binaryCoefficients = binaralizeCoefficientsViaAIC(coefficients, optimalGroupSize);
end