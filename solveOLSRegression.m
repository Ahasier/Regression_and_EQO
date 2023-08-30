function [estimatedCoefficients, trainingSets, testSets] = solveOLSRegression(abundanceData, functionalOutput, numPermutations, settings)
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

% Get the number of samples from the input abundance data
numberOfSamples = size(abundanceData, 1);

% Initialize variables to store intermediate results during permutations
currentBeta = zeros(size(abundanceData, 2), numPermutations);
trainingSets = zeros(numberOfSamples / 2, numPermutations);
testSets = zeros(numberOfSamples / 2, numPermutations);

% Loop over the specified number of permutations
for n = 1:numPermutations
    % Divide data into training and testing sets
    [trainingSets(:,n), testSets(:,n)] = generatePermutations(numberOfSamples);
    
    % Construct the optimization problem
    [regressionModel, x0] = formulateOptimization(abundanceData(trainingSets(:,n), :), functionalOutput(trainingSets(:,n)), [], 'OLS', settings);
    
    % Solve the optimization problem
    opts = optimoptions('ga');
    sol = solve(regressionModel, x0, 'Options', opts);
    
    % Get estimated betas on current training sample
    currentBeta(:, n) = sol.beta;
    
    % Handle extra features if needed
    % currentBeta(:, n) = handleExtraPhylogeneticFeatures(currentBeta(:, n), abundanceData, TaraNames, Idx, Ladd);
end
% Updata estimated coefficients
estimatedCoefficients = mean(currentBeta, 2);
end