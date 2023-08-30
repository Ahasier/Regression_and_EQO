function [estimatedCoefficients, trainingSets, testSets] = solveUnifiedRegression(abundanceData, functionalOutput, Lambda, numPermutations, regressionMethod, settings)
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
numberOfSamples = size(abundanceData, 1);
estimatedCoefficients = zeros(size(abundanceData, 2), length(Lambda));

for m = 1:length(Lambda)
    currentBeta = zeros(size(abundanceData, 2), numPermutations);
    trainingSets = zeros(numberOfSamples / 2, numPermutations);
    testSets = zeros(numberOfSamples / 2, numPermutations);
    
    for n = 1:numPermutations
        % Divide data into training and testing sets
        [trainingSets(:,n,m), testSets(:,n,m)] = generatePermutations(numberOfSamples);
        
        % Construct the optimization problem based on the current method and data split
        [regressionModel, x0] = formulateOptimization(abundanceData(trainingSets(:,n,m), :), functionalOutput(trainingSets(:,n,m)), Lambda(m), regressionMethod, settings);
        
        % Configure and solve the optimization problem using Genetic Algorithm (GA)
        opts = optimoptions('ga');
        sol = solve(regressionModel, x0, 'Options', opts);
        
        % Store the estimated coefficients for the current permutation
        currentBeta(:, n) = sol.beta;
        
        % (Optional) Adjust coefficients if considering extra phylogenetic features
        % currentBeta(:, n) = handleExtraPhylogeneticFeatures(currentBeta(:, n), TaraNames, Idx, Ladd);
    end
    
    % Average the estimated coefficients over all permutations for the current Lambda
    estimatedCoefficients(:, m) = mean(currentBeta, 2);
end
end