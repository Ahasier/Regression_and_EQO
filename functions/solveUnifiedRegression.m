function estimatedCoefficients = solveUnifiedRegression(trainingData, trainingOutput, regressionMethod, settings)
% SOLVEUNIFIEDREGRESSION solves regression problems using various methods
%
% INPUTS:
%   trainingData: Matrix representing the training data.
%   trainingOutput: Vector representing the training output.
%   regressionMethod: String specifying the regression method (e.g., 'LASSO', 'Ridge').
%   settings: A structure containing additional parameters and settings, including a range for the regularization parameter Lambda.
%
% OUTPUT:
%   estimatedCoefficients: Matrix of estimated regression coefficients for each value of Lambda.

% Initialize the size of the output matrix
Lambda = setLambdaRange(settings.maxLambda);
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
    if useExtraFeatures(settings)
        estimatedCoefficients(:, m) = handleExtraPhylogeneticFeatures(estimatedCoefficients(:, m), size(trainingData, 2), extraPhyloVars.Idx, extraPhyloVars.addedLeaves);
    end
end
end