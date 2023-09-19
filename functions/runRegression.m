function [binaryCoefficients, optimalGroupSize] = runRegression(trainingData, trainingOutput, regressionMethod, settings, varargin)
% Step (b) 1: Regression on training data
coefficients = computeRegression(trainingData, trainingOutput, regressionMethod, settings);

% If using regressions other than OLS (i.e, LASSO), cross-validate to get the optimal lambda
if ~strcmp(regressionMethod, 'OLS')
    testData = varargin{1};
    testOutput = varargin{2};
    [coefficients, ~] = crossValidationForLASSO(testData, testOutput, coefficients, settings.maxLambda);
end

% Step (b) 2: Depending on settings.Threshold, apply either AIC method, or
% a given threshold, or cross-validation over different threshold, to get
% the binary coefficients and optimal group size.
[groupingMethod, vararginToMethod] = determineGroupingMethod(settings);
[binaryCoefficients, optimalGroupSize] = groupingMethod(trainingData, trainingOutput, coefficients, vararginToMethod{:});
end

%% Helper functions
function [groupingMethod, vararginToMethod] = determineGroupingMethod(settings)
if isnan(settings.Threshold)
    groupingMethod = @runAIC;
    vararginToMethod = {};
else
    groupingMethod = @runThreshold;
    vararginToMethod = {settings};
end
end

function [binaryCoefficients, optimalGroupSize] = runThreshold(trainingData, trainingOutput, coefficients, settings)
% Use certain threshold method to get binary coefficients
[thresholdMethod, varargin] = determineThresholdMethod(settings);
binaryCoefficients = thresholdMethod(trainingData, trainingOutput, coefficients, varargin);

% Get group size from the binary coefficients
optimalGroupSize = getGroupSize(binaryCoefficients);
end

function [thresholdMethod, varargin] = determineThresholdMethod(settings)
if ischar(settings.Threshold) && strcmp(settings.Threshold, 'cv')
    thresholdMethod = @runThresholdCrossValidation;
    varargin = settings.Beta0;
elseif isnumeric(settings.Threshold)
    thresholdMethod = @runGivenThreshold;
    varargin = settings.Threshold;
end
end

function thresholedCoefficients = runGivenThreshold(~, ~, coefficients, threshold)
% Threshold coefficients according to the given threshold
thresholedCoefficients = coefficients > threshold;
end

function optimalGroupSize = getGroupSize(binaryCoefficients)
optimalGroupSize = sum(binaryCoefficients);
end