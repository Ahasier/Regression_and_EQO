function coefficients = runRegression(trainingData, trainingOutput, regressionMethod, settings, varargin)
% Step (b) 1: Regression on training data
coefficients = computeRegression(trainingData, trainingOutput, regressionMethod, settings);

% If using regressions other than OLS (i.e, LASSO), cross-validate to get the optimal lambda
if ~strcmp(regressionMethod, 'OLS')
    if length(varargin) == 2
        testData = varargin{1};
        testOutput = varargin{2};
    else
        testData = varargin{2};
        testOutput = varargin{3};
    end
    [coefficients, ~] = crossValidationForLASSO(testData, testOutput, coefficients, settings.maxLambda);
end

% Handle extra phylogenetic features if needed
if useExtraFeatures(settings)
    extraPhyloVars = varargin{1};
    coefficients = handleExtraPhylogeneticFeatures(coefficients, extraPhyloVars.numTaxa, extraPhyloVars.Idx, extraPhyloVars.addedLeaves);
end
end