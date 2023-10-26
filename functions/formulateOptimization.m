function [regressionModel, x0] = formulateOptimization(abundanceSubset, functionalOutputSubset, lambda, regressionMethod, settings)
% FORMULATEOPTIMIZATION formulates an optimization problem based on the specified regression
% method and returns a regression model (optimproblem object) and initial point x0
% 
% INPUTS:
%   abundanceSubset: Subset of the abundance data.
%   functionalOutputSubset: Subset of the functional output data.
%   lambda: Regularization parameter for the regression.
%   regressionMethod: String specifying the regression method (e.g., 'LASSO', 'Ridge', 'L0').
%   settings: A structure containing additional parameters and settings.
%
% OUTPUTS:
%   regressionModel: An optimization problem object suitable for the specified regression method.
%   x0: Initial point for the optimization problem.

% Detect is positivity requirement is needed
if isfield(settings, 'requirePositivity') && strcmp(settings.requirePositivity, 'On')
    requirePositivity = true;
else
    requirePositivity = false;
end

% Select the appropriate regression model based on the specified method
switch regressionMethod
    case 'LASSO'
        [regressionModel, x0] = lassoObjective(abundanceSubset, functionalOutputSubset, lambda, requirePositivity);
    case 'Ridge'
        [regressionModel, x0] = ridgeObjective(abundanceSubset, functionalOutputSubset, lambda, requirePositivity);
    case 'L0'
        [regressionModel, x0] = l0Objective(abundanceSubset, functionalOutputSubset, lambda, settings.RegPower);
    case 'Nonlinear'
        [regressionModel, x0] = nonlinearDesignedObjective(abundanceSubset, functionalOutputSubset, lambda, requirePositivity);
    case 'Combined'
        [regressionModel, x0] = combinedObjective(abundanceSubset, functionalOutputSubset, lambda, requirePositivity);
    case 'Binary'
        [regressionModel, x0] = binaryObjective(abundanceSubset, functionalOutputSubset);
    case 'OLS'
        [regressionModel, x0] = OLSObjective(abundanceSubset, functionalOutputSubset, requirePositivity);
end
end

%% Helper functions
% LASSO Regression Objective Function
% Formulate LASSO regression as an optimization problem
function [prob, x0] = lassoObjective(trainingData, functionalOutput, lambda, requirePositivity)
% Define optimization variable for LASSO coefficients
beta = defineOptimizationVariable(requirePositivity, trainingData);

% Compute the residual sum of squares
residual = trainingData * beta - functionalOutput;

% Define the LASSO objective function
obj = sum(residual.^2) + lambda * sum(beta.^2);

% Create the optimization problem
prob = optimproblem('Objective', obj);

% Define the initial point for the optimization algorithm
x0.beta = zeros(size(trainingData, 2), 1)';
end

% Ridge Regression Objective Function
% Formulate Ridge regression as an optimization problem
function [prob, x0] = ridgeObjective(trainingData, functionalOutput, lambda, requirePositivity)
% Define optimization variable for Ridge coefficients
beta = defineOptimizationVariable(requirePositivity, trainingData);

% Compute the residual sum of squares
residual = trainingData * beta - functionalOutput;

% Define the Ridge objective function
obj = sum(residual.^2) + lambda * sum(beta.^2);

% Create the optimization problem
prob = optimproblem('Objective', obj);

% Define the initial point for the optimization algorithm
x0.beta = zeros(size(trainingData, 2), 1)';
end

% L0 Regression Objective Function
% Formulate L0 norm regression as an optimization problem
function [prob, x0] = l0Objective(trainingData, functionalOutput, lambda, regPower)
% Define optimization variable for L0 regression coefficients
beta = optimvar('beta', size(trainingData, 2), 'LowerBound', 0);

% Compute the residual sum of squares
residual = trainingData * beta - functionalOutput;

% Define the L0 objective function with a power regularization term
obj = sum(residual.^2) + lambda * sum(beta.^regPower);

% Create the optimization problem
prob = optimproblem('Objective', obj);

% Define the initial point for the optimization algorithm
x0.beta = zeros(size(trainingData, 2), 1)';
end

% Binary Regression Objective Function
% Formulate binary regression as an optimization problem
function [prob, x0] = binaryObjective(trainingData, functionalOutput)
% Define optimization variable for binary regression coefficients
beta = optimvar('beta', size(trainingData, 2), 'Type', 'integer', 'LowerBound', 0, 'UpperBound', 1);

% Compute the residual sum of squares
residual = trainingData * beta - functionalOutput;

% Define the binary regression objective function
obj = sum(residual.^2);

% Create the optimization problem
prob = optimproblem('Objective', obj);

% Define the initial point for the optimization algorithm
x0.beta = zeros(size(trainingData, 2), 1)';
end

% Nonlinear Designed Regression Objective Function
% Formulate nonlinear designed regression as an optimization problem
function [prob, x0] = nonlinearDesignedObjective(trainingData, functionalOutput, lambda, requirePositivity)
threshold = 0.5;

% Define optimization variable for nonlinear designed regression coefficients
beta = defineOptimizationVariable(requirePositivity, trainingData);

% Compute the residual sum of squares
residual = trainingData * beta - functionalOutput;

% Define the nonlinear designed penalty term
obj = sum(residual.^2) + lambda*sum(nonlinearDesignedPenalty(beta, threshold));

% Create the optimization problem
prob = optimproblem('Objective', obj);

% Define the initial point for the optimization algorithm
x0.beta = zeros(size(trainingData, 2), 1)';
end

% Combined Regression Objective Function
% Formulate combined regression as an optimization problem
function [prob, x0] = combinedObjective(trainingData, functionalOutput, lambda, requirePositivity)
threshold = 0.5;

% Define optimization variable for L-combined regression coefficients
beta = defineOptimizationVariable(requirePositivity, trainingData);

% Compute the residual sum of squares
residual = trainingData * beta - functionalOutput;

% Define the L-combined penalty term
obj = sum(residual.^2) + sum(combinedPenalty(beta, lambda, threshold));

% Create the optimization problem
prob = optimproblem('Objective', obj);

% Define the initial point for the optimization algorithm
x0.beta = zeros(size(trainingData, 2), 1)';
end

% Ordinary Least Squares (OLS) Regression Objective Function
% Formulate OLS regression as an optimization problem
function [prob, x0] = OLSObjective(trainingData, functionalOutput, requirePositivity)
% Define optimization variable for OLS coefficients
beta = defineOptimizationVariable(requirePositivity, trainingData);

% Compute the residual sum of squares
residual = trainingData * beta - functionalOutput;

% Define the OLS objective function
obj = sum(residual.^2);

% Create the optimization problem
prob = optimproblem('Objective', obj);

% Define the initial point for the optimization algorithm
x0.beta = zeros(size(trainingData, 2), 1)';
end

function beta = defineOptimizationVariable(requirePositivity, trainingData)
if requirePositivity
    beta = optimvar('beta', size(trainingData, 2), 'LowerBound', 0);
else
    beta = optimvar('beta', size(trainingData, 2));
end
end

% Define the nonlinear designed penalty term for regression coefficients
function Penalty = nonlinearDesignedPenalty(beta, threshold)
% Compute the penalty based on the specified nonlinear function
Penalty = - log((exp((beta - threshold)/2)+exp(-(beta - threshold)/2))/2 - 0.99);
end

% Define the combined penalty term for regression coefficients
function Penalty = combinedPenalty(beta, lambda, threshold)
lambda2 = 3;
lambda3 = 3;
% Compute the penalty based on the specified combined function
Penalty = - log(((exp((beta - threshold).*lambda3) + exp(-(beta - threshold).*lambda3))/2 - 1).^lambda2.*exp(-lambda.*beta) + 0.01);
end