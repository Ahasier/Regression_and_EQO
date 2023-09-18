% Calculate the AIC values under different group sizes
function aicValues = evaluateAIC(trainingData, trainingOutput, coefficients)
% Initialization
numTaxa = length(coefficients);
aicValues = zeros(length(coefficients), 1);

% Sort coefficients by descend for later use
[~, sortedTaxaIndices] = sort(coefficients, 'descend');

% Loop over different group sizes and compute the corresponding AIC values
for n = 1:numTaxa
    aicValues(n) = computeAIC(n, numTaxa, trainingData, trainingOutput, sortedTaxaIndices);
end
end