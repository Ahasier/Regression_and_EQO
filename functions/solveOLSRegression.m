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

if isfield(settings, 'requirePositivity') && strcmp(settings.requirePositivity, 'On')
    % Construct the optimization problem
    [regressionModel, x0] = formulateOptimization(trainingData, trainingOutput, [], 'OLS', settings);
    
    % Solve the optimization problem
    opts = optimoptions('ga');
    sol = solve(regressionModel, x0, 'Options', opts);
    
    % Get estimated coefficients on current training samples
    estimatedCoefficients = sol.beta;
elseif size(trainingData, 1) == size(trainingData, 2)
    estimatedCoefficients = pinv(trainingData) * trainingOutput;
else
    estimatedCoefficients = trainingData \ trainingOutput;
end
end