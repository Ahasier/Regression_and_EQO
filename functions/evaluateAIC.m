% Calculate the AIC values under different group sizes
function aicc = evaluateAIC(trainingData, trainingOutput, coefficients)
% Initialization
% maxNumParameters = sum(coefficients > 0) + 1;
maxNumParameters = length(coefficients) + 1;
numSamples = length(trainingOutput);
lenAIC = min(maxNumParameters - 1, numSamples - 2);
aicc = zeros(lenAIC, 1);

% Sort coefficients by descend for later use
[~, sortedTaxaIndices] = sort(coefficients, 'descend');

% Loop over different group sizes and compute the corresponding AIC values
for n = 1:lenAIC
    % Get the capped coefficients
    cappedCoefficients = capCoefficients(coefficients, n, sortedTaxaIndices);
    
    % Compute the aic value of at this group size
    aicValue = computeAIC(n + 1, trainingData, trainingOutput, cappedCoefficients);
    
    % Compute the corrected AICc
    numParameters = n + 1;
    aicc(n) = aicValue + 2*(numParameters*(numParameters + 1))/(numSamples - numParameters - 1);
end
end

%% Helper function
function cappedCoefficients = capCoefficients(coefficients, n, sortedTaxaIndices)
idx = sortedTaxaIndices(n + 1:end);
% cappedCoefficients = coefficients;
cappedCoefficients = ones(size(coefficients));
cappedCoefficients(idx) = 0;
end