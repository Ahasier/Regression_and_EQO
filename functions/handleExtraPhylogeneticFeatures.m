% If phylogenetic information is incorporated, update estimated coefficients using extra features.
function updatedCoefficients = handleExtraPhylogeneticFeatures(coefficients, numTaxa, Idx, addedLeaves)
updatedAllBeta = coefficients(1:numTaxa);
for id = 1:length(Idx)
    updatedAllBeta(addedLeaves(id)) = updatedAllBeta(addedLeaves(id)) + coefficients(numTaxa + Idx(id));
end
updatedCoefficients = updatedAllBeta;
end