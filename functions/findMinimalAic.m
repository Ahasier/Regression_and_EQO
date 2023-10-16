% Find the optimal group size by minimizing AIC values
function [minimalAicValue, optimalGroupSize] = findMinimalAic(aicValues)
% Find the local minimal AICs
IndexCandidateMinAIC = findLocalMin(aicValues);

% Filter the candidate minimal AICs by ignoring those too close to the
% rightmost sharp drop region
% IndexCandidateMinAIC = filterMinAIC(aicValues, IndexCandidateMinAIC);

% Find the overal minimal AIC value
[minimalAicValue, optimalGroupSize] = findOveralMin(aicValues, IndexCandidateMinAIC);
end

%% Helper functions
function localMins = findLocalMin(data)
% Initialize an empty array to store the indices of local minimums
localMins = [];

% Check the length of the data
n = length(data);

% Iterate through the data starting from the second element and ending at the second last element
for i = 2:n-1
    % Check if the current element is less than its neighbors
    if data(i) <= data(i-1) && data(i) <= data(i+1)
        % If true, add the index to the localMins array
        localMins = [localMins, i];
    end
end
end

function filteredIndices = filterMinAIC(data, indices)
cutOff = length(data) - 20;
filteredIndices = indices(indices < cutOff);
end

function [globalMinValue, globalMinIndex] = findOveralMin(data, localMins)
% Find the global minimum from the local minimums
if isempty(localMins)
    [globalMinValue, globalMinIndex] = min(data);
    warning('No local minimal aic value detected because the aic curve is monotonous. Used the global minimum index at which %d taxa is selected.', globalMinIndex);
else
    [~, idx] = min(data(localMins));
    globalMinIndex = localMins(idx);
    globalMinValue = data(globalMinIndex);
end
end