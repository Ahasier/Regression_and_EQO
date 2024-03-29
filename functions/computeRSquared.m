% Compute the R2 values and Squared Error for given data and coefficients
function [R2, SSR] = computeRSquared(abundanceData, actualOutput, coefficients)
% Compute predicted values for the test set
predictedOutput = abundanceData * coefficients;
% Compute residuals
residuals = actualOutput - predictedOutput;
% % Centering the residuals
% residualsCentered = residuals - mean(residuals);
% Compute sum of squared residuals
SSR = sum((residuals).^2);
% Compute total sum of squares
TSS = sum((actualOutput - mean(actualOutput)).^2);
% Calculate R2 as 1 minus the ratio of SSR to TSS
R2 = 1 - (SSR / TSS);
end