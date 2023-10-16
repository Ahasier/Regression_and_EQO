function [binaryCoefficients, optimalGroupSize] = runAIC(trainingData, trainingOutput, coefficients)
% Step (b) 2: Calculate the AIC values under different group sizes to find the optimal one
aicValues = evaluateAIC(trainingData, trainingOutput, coefficients);
[~, optimalGroupSize] = findMinimalAic(aicValues);
% optimalGroupSize = detectMinimumBeforeDrop(aicValues, 3);

% Step (b) 3: Calculate the binary results sparsified from AIC step
binaryCoefficients = binaralizeCoefficientsViaAIC(coefficients, optimalGroupSize);
end