% Calculate the binary coefficients sparsified from AIC step
function binaryCoefficients = binaralizeCoefficientsViaAIC(coefficients, optimalGroupSize)
binaryCoefficients = ones(size(coefficients, 1), 1);

% Sort coefficients by descend for later use
[~, sortedTaxaIndices] = sort(coefficients, 'descend');

% Filter taxa by the optimal group size
IndicesInGroup = sortedTaxaIndices(1:optimalGroupSize);

% Set coefficients of other taxa to 0
for taxonIdx = 1:length(coefficients)
    if ~ismember(taxonIdx, IndicesInGroup)
        binaryCoefficients(taxonIdx) = 0;
    end
end
end