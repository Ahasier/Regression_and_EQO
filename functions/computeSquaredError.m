% Compute the R2 values and Squared Error for given data and coefficients
function SquaredError = computeSquaredError(abundanceData, actualOutput, coefficients)
% Compute predicted values for the test set
predictedOutput = abundanceData * coefficients;
% Compute residuals
residuals = actualOutput - predictedOutput;
% Compute mean squared error for the test set
SquaredError = sum(residuals.^2) / length(residuals);
end

