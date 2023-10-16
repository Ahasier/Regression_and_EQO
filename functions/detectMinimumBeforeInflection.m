function minIndex = detectMinimumBeforeInflection(data, windowSize)
% Step 1: Smooth the data using a moving average filter
smoothedData = movmean(data, windowSize);

% Step 2: Compute the second derivative of the smoothed data
secondDerivative = diff(smoothedData, 2);

% Step 3: Identify the inflection point
inflectionIndices = find(diff(sign(secondDerivative)) > 0);

if isempty(inflectionIndices)
    inflectionIndex = length(data) - 1; % -1 due to second derivative offset
else
    inflectionIndex = inflectionIndices(1);
end

% Step 4: Find the global minimum before the inflection point
[~, relativeMinIndex] = min(data(1:inflectionIndex));
minIndex = relativeMinIndex; % adjust for indexing from the start
end