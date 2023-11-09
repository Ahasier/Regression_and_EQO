function [groupedCoefficients, optimalAicValue, optimalGroupSize] = applyGroupingMethod(trainingData, trainingOutput, rawCoefficients, settings)
% Step (b) 2: Depending on settings.Threshold, apply either AIC method, or
% a given threshold, or cross-validation over different threshold, to get
% the binary coefficients and optimal group size.
if useExtraFeatures(settings)
    % Turning `trainingData` back to the ungrouped one for the next step
    trainingData = trainingData(:, 1:extraPhyloVars.numTaxa);
end
[groupingMethod, vararginToMethod] = determineGroupingMethod(settings);
[groupedCoefficients, optimalAicValue, optimalGroupSize] = groupingMethod(trainingData, trainingOutput, rawCoefficients, vararginToMethod{:});
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

function [binaryCoefficients, optimalAicValue, optimalGroupSize] = runThreshold(trainingData, trainingOutput, coefficients, settings)
% Use certain threshold method to get binary coefficients
[thresholdMethod, varargin] = determineThresholdMethod(settings);
binaryCoefficients = thresholdMethod(trainingData, trainingOutput, coefficients, varargin);

% Get group size from the binary coefficients
optimalGroupSize = getGroupSize(binaryCoefficients);
optimalAicValue = computeAIC(optimalGroupSize, trainingData, trainingOutput, coefficients(binaryCoefficients));
end

function [thresholdMethod, varargin] = determineThresholdMethod(settings)
if ischar(settings.Threshold) && strcmp(settings.Threshold, 'cv')
    thresholdMethod = @runThresholdCrossValidation;
    varargin = settings.Beta0;
elseif isnumeric(settings.Threshold)
    thresholdMethod = @runGivenThreshold;
    varargin = settings.Threshold;
else
    error('Threshold method must either be "cv" or a given numeric value, or nan.');
end
end

function thresholedCoefficients = runGivenThreshold(~, ~, coefficients, threshold)
% Threshold coefficients according to the given threshold
thresholedCoefficients = coefficients > threshold;
end

function optimalGroupSize = getGroupSize(binaryCoefficients)
optimalGroupSize = sum(binaryCoefficients);
end