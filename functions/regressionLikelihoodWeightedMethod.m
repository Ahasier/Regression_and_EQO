function aicValue = regressionLikelihoodWeightedMethod(n, trainingData, trainingOutput, coefficients, settings)
if isfield(settings, 'weight') && strcmp(settings.weight, 'On')
%     % Find optimal weight by minimizing loglikelihood
%     weight = findOptimalWeight(trainingOutput, trainingData, coefficients);
%     
%     % Calculate the minimal loglikelihoood
%     coefficients = weight.*coefficients;
%     y_pred = trainingData * coefficients;
%     residuals = trainingOutput - y_pred;
%     logLikelihood = regressionLikelihood(residuals);
%     
%     % Compute AIC value
%     aicValue = 2*n - 2*logLikelihood;
    
    % Loop over different weights to pick an optimal one
    weights = setWeightsRange();
    lenWeights = length(weights);
    logLikelihoods = zeros(1, lenWeights);
    
    for w = 1:lenWeights
        weight = weights(w);
        
        % Weight the coefficients
        weightedCoefficients = weight.*coefficients;
        
        % Compute predicted values
        y_pred = trainingData * weightedCoefficients;
        
        % Compute residuals
        residuals = trainingOutput - y_pred;
        
        % Compute Log-likelihood
        logLikelihoods(w) = regressionLikelihood(residuals);
    end
    
    % Find the optimal weight that maximizes Log-likelihood
    maxLogLikelihood = max(logLikelihoods);
    
    % Compute AIC value
    aicValue = 2*n - 2*maxLogLikelihood;
else
    % Compute predicted values
    y_pred = trainingData * coefficients;
    
    % Compute residuals
    residuals = trainingOutput - y_pred;
    
    % Compute Log-likelihood
    logLikelihood = regressionLikelihood(residuals);
    
    % Compute AIC value
    aicValue = 2*n - 2*logLikelihood;
end
end

function w = findOptimalWeight(y, X, beta)
n = length(y);
y_pred = X * beta;
avr_y = mean(y);
avr_y_pred = mean(y_pred);
w = (y' * y_pred - (n/2) * avr_y * avr_y_pred)/(y_pred' * y_pred - (n/2) * avr_y_pred.^2);
end