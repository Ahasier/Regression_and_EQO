function [optimalCoefficients, resultsForDiagnostics] = computeRegressionAndCrossValidationViaAIC(abundanceData, functionalOutput, numPermutations, regressionMethod, settings)
% Initialize variables to store AIC values and corresponding coefficients
[allBinaryCoefficients, allAicValues, allOutSampleR2] = initializeAICCrossValidationResults();
if ~strcmp(regressionMethod, 'OLS')
    optimalLambda = zeros(1,numPermutations);
end

for i = 1:numPermutations
    % Step (a): Split data
    [trainingData, testData, trainingOutput, testOutput] = splitData(abundanceData, functionalOutput);

    % Step (b) 1: Regression on training data
    coefficients = computeRegression(trainingData, trainingOutput, regressionMethod, settings);

    % If using regressions other than OLS (i.e, LASSO), cross-validate to get the optimal lambda
    if ~strcmp(regressionMethod, 'OLS')
        [coefficients, optimalLambda(i)] = crossValidationForLASSO(testData, testOutput, coefficients, settings.maxLambda);
    end
    
    % Step (b) 2: Calculate the AIC values under different group sizes to find the optimal one
    aicValues = evaluateAIC(trainingData, trainingOutput, coefficients);
    [minimalAicValue, optimalGroupSize] = findMinimalAic(aicValues);
    
    % Step (c): Calculate the binary results (sparsified from AIC step) and 
    % cross-validation R^2 on test subset
    binaryCoefficients = binaralizeCoefficientsViaAIC(coefficients, optimalGroupSize);
    R2OutSamples = computeRSquared(testData, testOutput, binaryCoefficients);
    
    % Store results
    allBinaryCoefficients(:, i) = binaryCoefficients;
    allOutSampleR2(i) = R2OutSamples;
    allAicValues(i,:) = aicValues;
end

% Step (d) 1: Calculate the cumulative R^2 for each taxon
cumulativeR2 = calculateCumulativeR2(numPermutations, abundanceData, allBinaryCoefficients, allOutSampleR2);

% Step (d) 2: Get the best coefficients by selecting the top k = optimalGroupSize taxa based on cumulative R^2
optimalCoefficients = selectTopKTaxa(cumulativeR2, optimalGroupSize, allBinaryCoefficients);

% Store results as a structure
if ~strcmp(regressionMethod, 'OLS')
    resultsForDiagnostics = struct('allBinaryCoefficients', allBinaryCoefficients, 'allOutSampleR2', allOutSampleR2, 'minimalAicValue', minimalAicValue, 'cumulativeR2', cumulativeR2);
else
    resultsForDiagnostics = struct('allBinaryCoefficients', allBinaryCoefficients, 'allOutSampleR2', allOutSampleR2, 'optimalLambda', optimalLambda, 'minimalAicValue', minimalAicValue, 'cumulativeR2', cumulativeR2);
end
end

%% Helper functions
function [allCoefficients, aicValues, allOutSampleR2] = initializeAICCrossValidationResults(abundanceData, numPermutations)
allCoefficients = zeros(size(abundanceData, 2), numPermutations);
allOutSampleR2 = zeros(1, numPermutations);
aicValues = zeros(size(abundanceData, 2), numPermutations);
end

% Calculate the AIC values under different group sizes
function aicValue = evaluateAIC(trainingData, trainingOutput, coefficients)
% Initialization
aicValue = zeros(1,length(coefficients));

% Sort coefficients by descend for later use
[~, sortedTaxaIndices] = sort(coefficients, 'descend');

% Loop over different group sizes and compute the corresponding AIC values
for n = 1:length(coefficients)
    idx = sortedTaxaIndices(1:n);
    groupAssemblage = zeros(length(coefficients), 1);
    groupAssemblage(idx) = 1;
    
    [r2, ssr] = computeRSquared(trainingData, trainingOutput, groupAssemblage);
    
    k = size(trainingData, 2) + 1; % Number of predictors + intercept
    L = exp(-ssr/2);
    aic = 2*k - 2*log(L);
    
    % Store the AIC value
    aicValue(n) = aic;
end
end

% Find the optimal group size by minimizing AIC values
function [minimalAicValue, optimalGroupSize] = findMinimalAic(aicValues)
avgAicValues = mean(aicValues);
[minimalAicValue, optimalGroupSize] = min(avgAicValues);
end

% Calculate the binary coefficients sparsified from AIC step
function binaryCoefficients = binaralizeCoefficientsViaAIC(coefficients, optimalGroupSize)
binaryCoefficients = ones(size(allCoefficients, 1), 1);

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

% Step 6: Calculate the cumulative R^2 for each taxon
function cumulativeR2 = calculateCumulativeR2(numPermutations, abundanceData, allCoefficients, allOutSampleR2)
cumulativeR2 = zeros(size(abundanceData, 2), 1);
for taxonIdx = 1:size(abundanceData, 2)
    for permIdx = 1:numPermutations
        cumulativeR2(taxonIdx) = cumulativeR2(taxonIdx) + allCoefficients(taxonIdx, permIdx) * allOutSampleR2(permIdx);
    end
end
end

% Get bet best coefficients by selecting the top k taxa based on cumulative R^2
function optimalCoefficients = selectTopKTaxa(cumulativeR2, k, allCoefficients)
optimalCoefficients = ones(size(allCoefficients, 1), 1);

% Select the top k taxa based on cumulative R^2
[~, sortedTaxaIndices] = sort(cumulativeR2, 'descend');
topKTaxaIndices = sortedTaxaIndices(1:k);

% Set coefficients of other taxa to 0
for taxonIdx = 1:size(allCoefficients, 1)
    if ~ismember(taxonIdx, topKTaxaIndices)
        optimalCoefficients(taxonIdx) = 0;
    end
end
end