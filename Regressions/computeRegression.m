function estimatedCoefficients = computeRegression(trainingData, trainingOutput, regressionMethod, settings)
% COMPUTEREGRESSIONS Performs regression based on input parameters.
% INPUTS:
%   trainingData: Matrix containing the training data.
%   trainingOutput: Vector containing the corresponding training output.
%   regressionMethod: String specifying the regression method ('L0', 'OLS', etc.).
%   settings: A structure with additional parameters and settings.
% 
% OUTPUT:
%   estimatedCoefficients: Vector of estimated regression coefficients based on the chosen regression method.

% Solve the regression problem based on the method provided
if strcmp(regressionMethod, 'L0')
    estimatedCoefficients = handleL0Regression(trainingData, trainingOutput, settings);
elseif strcmp(regressionMethod, 'OLS')
    estimatedCoefficients = solveOLSRegression(trainingData, trainingOutput, settings);
else
    estimatedCoefficients = solveUnifiedRegression(trainingData, trainingOutput, regressionMethod, settings);
end
end

%% Helper functions
function estimatedCoefficients = handleL0Regression(trainingData, trainingOutput, settings)
Lambda = 1:1:settings.maxLambda;
if isfield(settings, 'L0Option')
    if strcmp(settings.L0Option, 'MIQP')
        estimatedCoefficients = solveL0RegressionMIQP(trainingData, trainingOutput, Lambda);
    elseif strcmp(settings.L0Option, 'IHT')
        estimatedCoefficients = solveL0RegressionIterativeHardThresholding(trainingData, trainingOutput, selectedNames, Lambda, settings.T0, settings.alpha, settings.max_iter);
    end
end
end