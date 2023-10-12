function aicValue = computeAIC(n, trainingData, trainingOutput, coefficients)
% Compute the residuals
residuals = residualsOfPrediction(coefficients, trainingData, trainingOutput);

% Calculate the log-likelihood
logLikelihood = regressionLikelihood(residuals);

% Calculate the aic value
aicValue = 2*n - 2*logLikelihood;
end

%% Helper functions
function residuals = residualsOfPrediction(coefficients, abundance, y)
% P is the unweighted prediction of the outputs
P = abundance * coefficients;

% The weighted predicted outpus are `y_pred = a + b * P`, where the values 
% of a and b are inferred from a regression.
% A = [ones(length(P), 1), P];
A = P;
z = A \ y;
y_pred = A * z;

% Calculate the residual values
residuals = y - y_pred;
end