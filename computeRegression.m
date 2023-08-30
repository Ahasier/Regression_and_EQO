function [estimatedCoefficients, trainingSets, testSets] = computeRegression(abundanceData, functionalOutput, numPermutations, regressionMethod, settings)
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
    [estimatedCoefficients, trainingSets, testSets] = handleL0Regression(abundanceData, functionalOutput, numPermutations, settings);
elseif strcmp(regressionMethod, 'OLS')
    [estimatedCoefficients, trainingSets, testSets] = handleOLSRegression(abundanceData, functionalOutput, numPermutations, settings);
else
    [estimatedCoefficients, trainingSets, testSets] = handleUnifiedRegression(abundanceData, functionalOutput, numPermutations, regressionMethod, settings);
end
end

%% Helper functions
function [estimatedCoefficients, trainingSets, testSets] = handleL0Regression(abundanceData, functionalOutput, numPermutations, settings)
Lambda = 1:1:settings.maxLambda;
if isfield(settings, 'L0Option')
    if strcmp(settings.L0Option, 'MIQP')
        [estimatedCoefficients, trainingSets, testSets] = solveL0RegressionMIQP(abundanceData, functionalOutput, Lambda, numPermutations);
    elseif strcmp(settings.L0Option, 'IHT')
        [estimatedCoefficients, trainingSets, testSets] = solveL0RegressionIterativeHardThresholding(abundanceData, functionalOutput, selectedNames, Lambda, numPermutations, settings.T0, settings.alpha, settings.max_iter);
    end
end
end

function [estimatedCoefficients, trainingSets, testSets] = handleOLSRegression(abundanceData, functionalOutput, numPermutations, settings)
[estimatedCoefficients, trainingSets, testSets] = solveOLSRegression(abundanceData, functionalOutput, numPermutations, settings);
end

function [estimatedCoefficients, trainingSets, testSets] = handleUnifiedRegression(abundanceData, functionalOutput, numPermutations, regressionMethod, settings)
Lambda = 0:0.1:settings.maxLambda;
[estimatedCoefficients, trainingSets, testSets] = solveUnifiedRegression(abundanceData, functionalOutput, Lambda, numPermutations, regressionMethod, settings, Lambda);
end