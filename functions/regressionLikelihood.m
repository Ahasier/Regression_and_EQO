function logLikelihood = regressionLikelihood(residuals)
% X: matrix of predictors (each row is an observation, each column is a predictor)
% y: column vector of actual outputs
% beta: column vector of coefficients

% % Compute predicted values
% y_pred = X * beta;
% 
% % Compute residuals
% residuals = y - y_pred;

% Estimate variance of residuals
sigma2 = var(residuals);

% Compute log likelihood for each observation
logLikelihood_i = -0.5 * log(2 * pi * sigma2) - ((residuals.^2) / (2 * sigma2));

% Sum to get overall log likelihood
logLikelihood = sum(logLikelihood_i);
end
