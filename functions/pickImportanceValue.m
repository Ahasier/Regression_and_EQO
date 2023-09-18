function optimalCoefficients = pickImportanceValue(importanceValues)
numTaxa = length(importanceValues);
optimalCoefficients = zeros(numTaxa, 1);

% Select the the taxa with their importanceValue > 0.5, otherwise unselect
optimalCoefficients(importanceValues > 0.5) = 1;
end