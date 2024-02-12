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

% Include interception term in the optimization model
abundanceWithIntercept = [trainingData, ones(size(trainingData, 1), 1)];

if isfield(settings, 'requirePositivity') && strcmp(settings.requirePositivity, 'On')
    % Construct the optimization problem
    [regressionModel, x0] = formulateOptimization(abundanceWithIntercept, trainingOutput, [], 'OLS', settings);
    
    % Solve the optimization problem
    opts = optimoptions('ga');
    sol = solve(regressionModel, x0, 'Options', opts);
    
    % Get estimated coefficients on current training samples
    estimatedCoefficients = sol.beta;
elseif size(abundanceWithIntercept, 1) == size(abundanceWithIntercept, 2)
    estimatedCoefficients = pinv(abundanceWithIntercept) * trainingOutput;
else
    estimatedCoefficients = abundanceWithIntercept \ trainingOutput;
end

estimatedCoefficients = estimatedCoefficients(1:end - 1);
end