function estimatedCoefficients = solveOLSRegression(trainingData, trainingOutput, settings)
% SOLVEOLSREGRESSION performs multiple permutations of OLS regression to estimate the 
% association between microbial abundance data and a functional output. 
%
% Inputs:
%   abundanceData: matrix containing taxa abundance data (samples x taxa)
%   functionalOutput: vector containing the functional outputs (samples x 1)
%   numPermutations: the number of times the regression is to be permuted
%   regressionMethod: a string that specifies the type of regression method to be used
%   settings: a struct containing various settings for the optimization problem
%
% Outputs:
%   estimatedCoefficients: vector containing the estimated coefficients after permutations (taxa x 1)
%   trainingSets: indices of training samples for each permutation (training samples x numPermutations)
%   testSets: indices of testing samples for each permutation (testing samples x numPermutations)

% Construct the optimization problem
[regressionModel, x0] = formulateOptimization(trainingData, trainingOutput, [], 'OLS', settings);

% Solve the optimization problem
opts = optimoptions('ga');
sol = solve(regressionModel, x0, 'Options', opts);

% Get estimated coefficients on current training samples
estimatedCoefficients = sol.beta;

% Handle extra features if needed
% estimatedCoefficients = handleExtraPhylogeneticFeatures(estimatedCoefficients, abundanceData, TaraNames, Idx, Ladd);
end