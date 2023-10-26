function logLikelihood = regressionLikelihood(residuals)
% X: matrix of predictors (each row is an observation, each column is a predictor)
% y: column vector of actual outputs
% beta: column vector of coefficients

% % Compute predicted values
% y_pred = X * beta;
% 
% % Compute residuals
% residuals = y - y_pred;

% n is the number of observables
n = length(residuals);

% Estimate variance of residuals
% RSS = sum(residuals.^2);
% sigma2 = RSS/n;
sigma2 = rms(residuals)^2;
% sigma2 = var(residuals)*(n - 1)/n;

% Compute log likelihood
logLikelihood = - n/2 * (log(2 * pi * sigma2) + 1);
% logLikelihood = - n/2 * log(2 * pi * sigma^2) - sum(((residuals.^2) / (2 * sigma^2)));
end
