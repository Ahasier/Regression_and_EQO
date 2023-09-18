% Find the optimal group size by minimizing AIC values
function [minimalAicValue, optimalGroupSize] = findMinimalAic(aicValues)
[minimalAicValue, optimalGroupSize] = min(aicValues);
end