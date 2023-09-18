% Get bet best coefficients by selecting the top k taxa based on cumulative R^2
function optimalCoefficients = selectTopKTaxa(cumulativeR2, k)
numTaxa = length(cumulativeR2);
optimalCoefficients = ones(numTaxa, 1);

% Select the top k taxa based on cumulative R^2
[~, sortedTaxaIndices] = sort(cumulativeR2, 'descend');
topKTaxaIndices = sortedTaxaIndices(1:k);

% Set coefficients of other taxa to 0
for taxonIdx = 1:numTaxa
    if ~ismember(taxonIdx, topKTaxaIndices)
        optimalCoefficients(taxonIdx) = 0;
    end
end
end