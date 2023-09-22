function [binaryCoefficients, optimalGroupSize] = runAIC(trainingData, trainingOutput, coefficients, settings)
% Step (b) 2: Calculate the AIC values under different group sizes to find the optimal one
aicValues = evaluateAIC(trainingData, trainingOutput, coefficients, settings);
[~, optimalGroupSize] = findMinimalAic(aicValues);

% Step (b) 3: Calculate the binary results sparsified from AIC step
binaryCoefficients = binaralizeCoefficientsViaAIC(coefficients, optimalGroupSize);
end