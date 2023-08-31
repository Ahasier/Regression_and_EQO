function estimatedCoefficients = solveUnifiedRegression(trainingData, trainingOutput, regressionMethod, settings)
% SOLVEUNIFIEDREGRESSION Solves regression problems using various methods, with data splitting and averaging.
%
% INPUTS:
%   abundanceData: Matrix representing the abundance data.
%   functionalOutput: Vector representing the functional output of the system.
%   Lambda: Vector of regularization parameters.
%   numPermutations: Number of times to randomly permute the data for cross-validation.
%   regressionMethod: String specifying the regression method to be used ('LASSO', 'Ridge', etc.).
%   settings: A structure containing additional parameters and settings.
%
% OUTPUTS:
%   estimatedCoefficients: Matrix of estimated regression coefficients.
%   testSets: Matrix containing the test sets used in each permutation.

% Initialize the size of the output matrix
Lambda = 0:0.1:settings.maxLambda;
estimatedCoefficients = zeros(size(trainingData, 2), length(Lambda));

for m = 1:length(Lambda)
    % Construct the optimization problem based on the current method and data split
    [regressionModel, x0] = formulateOptimization(trainingData, trainingOutput, Lambda(m), regressionMethod, settings);
    
    % Configure and solve the optimization problem using Genetic Algorithm (GA)
    opts = optimoptions('ga');
    sol = solve(regressionModel, x0, 'Options', opts);
    
    % Store the estimated coefficients for the current Lambda on current training samples
    estimatedCoefficients(:, m) = sol.beta;
    
    % (Optional) Adjust coefficients if considering extra phylogenetic features
    % estimatedCoefficients(:, m) = handleExtraPhylogeneticFeatures(estimatedCoefficients(:, m), TaraNames, Idx, Ladd);
end
end