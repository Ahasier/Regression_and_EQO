function estimatedCoefficients = solveOLSRegression(trainingData, trainingOutput, settings)
% SOLVEOLSREGRESSION performs OLS regression to estimate the 
% association between microbial abundance data and a functional output. 
%
% INPUTS:
%   trainingData: Matrix containing the training data.
%   trainingOutput: Vector representing the training output.
%   settings: A structure containing various settings for the optimization problem.
%
% OUTPUT:
%   estimatedCoefficients: Vector containing the estimated regression coefficients after permutations.

if strcmp(settings.requirePositivity, 'Off')
    estimatedCoefficients = trainingData \ trainingOutput;
else
    % Construct the optimization problem
    [regressionModel, x0] = formulateOptimization(trainingData, trainingOutput, [], 'OLS', settings);
    
    % Solve the optimization problem
    opts = optimoptions('ga');
    sol = solve(regressionModel, x0, 'Options', opts);
    
    % Get estimated coefficients on current training samples
    estimatedCoefficients = sol.beta;
end
end