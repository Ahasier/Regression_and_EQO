function aicValue = computeAIC(k, trainingData, trainingOutput, coefficients)
% Compute the residuals
residuals = residualsOfPrediction(coefficients, trainingData, trainingOutput);

% Calculate the log-likelihood
% logLikelihood = regressionLikelihood(residuals);

% Calculate the aic value
% aicValue = 2*k - 2*logLikelihood;
n = length(residuals);
RSS = sum(residuals.^2);
aicValue = 2*k + n*log(RSS);
end

%% Helper functions
function residuals = residualsOfPredictionRescaleAll(coefficients, abundance, y)
% P is the unweighted prediction of the outputs
P = abundance * coefficients;

% The weighted predicted outpus are `y_pred = a + b * P`, where the values 
% of a and b are inferred from a regression.
A = [ones(length(P), 1), P];
% A = P;
z = A \ y;
y_pred = A * z;

% Calculate the residual values
residuals = y - y_pred;
end

function residuals = residualsOfPrediction(coefficients, abundance, y)
% Reduce the parameters of the model to only count for allowed taxa
c = coefficients ~= 0;
a = abundance(:,c);

A = [ones(size(a,1), 1), a];

% The coefficients of the reduced model are recovered from regression
b = A \ y;

% calculate y_predicted
y_pred = A * b;

% Calculate the residual values
residuals = y - y_pred;
end