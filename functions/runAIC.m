function [binaryCoefficients, optimalAicValue, optimalGroupSize] = runAIC(trainingData, trainingOutput, coefficients)
% Step (b) 2: Calculate the AIC values under different group sizes to find the optimal one
aicValues = evaluateAIC(trainingData, trainingOutput, coefficients);
[optimalAicValue, optimalGroupSize] = findMinimalAic(aicValues);
% optimalGroupSize = detectMinimumBeforeDrop(aicValues, 3);

% Step (b) 3: Calculate the binary results sparsified from AIC step
binaryCoefficients = binaralizeCoefficientsViaAIC(coefficients, optimalGroupSize);
end