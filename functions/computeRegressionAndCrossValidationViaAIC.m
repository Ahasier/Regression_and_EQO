function [optimalCoefficients, resultsForDiagnostics] = computeRegressionAndCrossValidationViaAIC(abundanceData, functionalOutput, numPermutations, regressionMethod, settings)
% Initialize variables to store AIC values and corresponding coefficients
[allCoefficients, aicValues, allOutSampleR2] = initializeAICCrossValidationResults();
if ~strcmp(regressionMethod, 'OLS')
    optimalLambda = zeros(1,numPermutations);
end

for i = 1:numPermutations
    % Step 1: Split data
    [trainingData, testData, trainingOutput, testOutput] = splitData(abundanceData, functionalOutput);
    
    % Step 2: Regression on training data
    coefficients = computeRegression(trainingData, trainingOutput, regressionMethod, settings);
    
    % If using regressions other than OLS (i.e, LASSO), cross-validate to get the optimal lambda
    if ~strcmp(regressionMethod, 'OLS')
        [coefficients, optimalLambda(i)] = crossValidationForLASSO(testData, testOutput, coefficients, settings.maxLambda);
    end
    
    % Step 3: Calculate the AIC values under different group sizes
    aicValue = evaluateAIC(trainingData, trainingOutput, coefficients);
    
    % Step 4: Calculate the cross-validation R^2 on test subset
    R2OutSamples = computeRSquared(testData, testOutput, coefficients);
    
    % Step 5: Store results
    allCoefficients(:, i) = coefficients;
    allOutSampleR2(i) = R2OutSamples;
    aicValues(i,:) = aicValue;
end

% Find the optimal group size k = optimalGroupSize by minimizing AIC
[minimalAicValue, optimalGroupSize] = findMinimalAIC(aicValues);

% Step 6: Calculate the cumulative R^2 for each taxon
cumulativeR2 = calculateCumulativeR2(numPermutations, abundanceData, allCoefficients, allOutSampleR2);

% Get bet best coefficients by selecting the top k = optimalGroupSize taxa based on cumulative R^2
optimalCoefficients = crossValidationViaAIC(cumulativeR2, optimalGroupSize, allCoefficients, abundanceData);

% Store results as a structure
if ~strcmp(regressionMethod, 'OLS')
    resultsForDiagnostics = struct('optimalLambda', optimalLambda, 'allOutSampleR2', allOutSampleR2, 'cumulativeR2', cumulativeR2);
else
    resultsForDiagnostics = struct('allOutSampleR2', allOutSampleR2, 'minimalAicValue', minimalAicValue, 'cumulativeR2', cumulativeR2);
end
end

%% Helper functions
function [allCoefficients, aicValues, allOutSampleR2] = initializeAICCrossValidationResults(abundanceData, numPermutations)
allCoefficients = zeros(size(abundanceData, 2), numPermutations);
allOutSampleR2 = zeros(1, numPermutations);
aicValues = zeros(size(abundanceData, 2), numPermutations);
end

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

function cumulativeR2 = calculateCumulativeR2(numPermutations, abundanceData, allCoefficients, allOutSampleR2)
cumulativeR2 = zeros(size(abundanceData, 2), 1);
for taxonIdx = 1:size(abundanceData, 2)
    for permIdx = 1:numPermutations
        cumulativeR2(taxonIdx) = cumulativeR2(taxonIdx) + allCoefficients(taxonIdx, permIdx) * allOutSampleR2(permIdx);
    end
end
end

function [minimalAicValue, optimalGroupSize] = findMinimalAIC(aicValues)
avgAicValues = mean(aicValues);
[minimalAicValue, optimalGroupSize] = min(avgAicValues);
end

% Get bet best coefficients by selecting the top k taxa based on cumulative R^2
function optimalCoefficients = crossValidationViaAIC(cumulativeR2, k, allCoefficients, abundanceData)
optimalCoefficients = zeros(size(allCoefficients, 1), 1);

% Select the top k taxa based on cumulative R^2
[~, sortedTaxaIndices] = sort(cumulativeR2, 'descend');
topKTaxaIndices = sortedTaxaIndices(1:k);

% Set coefficients of other taxa to 0
for taxonIdx = 1:size(abundanceData, 2)
    if ~ismember(taxonIdx, topKTaxaIndices)
        optimalCoefficients(taxonIdx) = 0;
    end
end
end