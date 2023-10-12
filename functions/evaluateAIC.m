% Calculate the AIC values under different group sizes
function aicValues = evaluateAIC(trainingData, trainingOutput, coefficients)
% Initialization
numTaxa = length(coefficients);
aicValues = zeros(length(coefficients), 1);

% Sort coefficients by descend for later use
[~, sortedTaxaIndices] = sort(coefficients, 'descend');

% Loop over different group sizes and compute the corresponding AIC values
for n = 1:numTaxa
    % Get the capped coefficients
    cappedCoefficients = capCoefficients(coefficients, n, sortedTaxaIndices);
    
    % Compute the aic value of at this group size
    aicValues(n) = computeAIC(n, trainingData, trainingOutput, cappedCoefficients);
%     aicValues(n) = regressionLikelihoodWeightedMethod(n, trainingData, trainingOutput, cappedCoefficients, settings);
end
end

%% Helper function
function cappedCoefficients = capCoefficients(coefficients, n, sortedTaxaIndices)
idx = sortedTaxaIndices(n + 1:end);
cappedCoefficients = coefficients;
cappedCoefficients(idx) = 0;
end