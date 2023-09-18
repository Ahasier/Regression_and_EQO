function bestCoefficients = runThresholdCrossValidation(trainingData, trainingOutput, coefficients, beta0)
% 1. Set threshold values for cross validation to go over
crossValidateThresholds = setCrossValidationThreshold(beta0);

% 2. Cross validate over all possible thresholds and get all thresholded coefficients
[allThresholdedCoefficients, allOutSampleErrors, ~] = handleCrossValidationForOLS(crossValidateThresholds,trainingData, trainingOutput, coefficients);

% 3. Find the threshold that gives the minimum out-of-sample error
[~, idx] = min(allOutSampleErrors(:));
bestCoefficients = allThresholdedCoefficients(:, idx);
end