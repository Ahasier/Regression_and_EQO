function estimatedCoefficients = computeRegression(trainingData, trainingOutput, regressionMethod, settings)
% COMPUTEREGRESSIONS Performs regression based on input parameters.
% INPUTS:
%   numPermutations: Number of times to randomly permute the data for cross-validation.
%   regressionMethod: String specifying the regression method to be used ('LASSO', 'Ridge', etc.).
%   settings: A structure containing additional parameters and settings.
% 
% OUTPUTS:
%   estimatedCoefficients: Matrix of estimated regression coefficients.
%   testSets: Matrix containing the test sets used in each permutation.

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