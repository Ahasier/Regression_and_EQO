function [TCM, importanceValues, crossValidatedCoefficients] = computeRegressionAndCrossValidation(abundanceData, functionalOutput, numPermutations, regressionMethod, settings)
% COMPUTEREGRESSIONANDCROSSVALIDATION integrates the regression and cross-validation processes.
% 
% INPUTS:
%   abundanceData: Matrix containing the abundance data of different taxa across samples.
%   functionalOutput: Vector containing the specific functional output in the samples.
%   numPermutations: Number of times to perform the integrated process.
%   regressionMethod: String specifying the regression method ('L0', 'OLS', etc.).
%   settings: A structure with additional parameters and settings.
%
% OUTPUTS:
%   Various metrics related to the integrated regression and cross-validation process, including coefficients, out-of-sample errors, R-squared values, and optimal thresholds.

% 1. Get table of candidate models (TCM)
[TCM, optimalGroupSize] = getTCM(abundanceData, functionalOutput, numPermutations, regressionMethod, settings);

% Check if threshold method is 'cv'
if ischar(settings.Threshold) && strcmp(settings.Threshold, 'cv')
    % If already used optimal threshold through cross-validation, no longer
    % need to go through extra crpss-validation step.
    importanceValues = [];
    crossValidatedCoefficien  ts = median(TCM{1:end-1, :}, 2);
else
    % 2. Get importance values from TCM
    importanceValues = getImportanceValues(TCM);
    
    % 3. Get the cross-validated coefficients
    % crossValidatedCoefficients = pickImportanceValue(importanceValues);
    crossValidatedCoefficients = selectTopKTaxa(importanceValues, optimalGroupSize);
end
end

%% Helper functions
function [TCM, optimalGroupSize] = getTCM(abundanceData, functionalOutput, numPermutations, regressionMethod, settings)
% Initialize TCM and other variables
numTaxa = size(abundanceData, 2);
TCM = initializeTCM(numTaxa, numPermutations);
optimalGroupSizes = zeros(1, numPermutations);

% Determine which method to use based on settings.Threshold
[trainingMethod, varargin, addTestData] = determineMethod(regressionMethod, settings);

for i = 1:numPermutations
    % Step (a): Split data
    [trainingData, testData, trainingOutput, testOutput] = splitData(abundanceData, functionalOutput);
    
    % Check if testData and testOutput need to be added to varargin
    if addTestData
        [varargin{end + 1:end + 2}] = deal(testData, testOutput);
    end
    
    % Step (b): Depending on the threshold value, use regressions or EQO to train on training data
    [coefficients, optimalGroupSizes(i)] = trainingMethod(trainingData, trainingOutput, varargin{:});
    
    % Step (c): Calculate the out of sample R^2 on test subset
    R2OutSamples = computeRSquared(testData, testOutput, coefficients);
    
    % Store results
    TCM{:, i} = [coefficients; R2OutSamples];
end

% Determine the optimal group size
optimalGroupSize = mode(optimalGroupSizes);
end

function TCM = initializeTCM(numTaxa, numPermutations)
% Initialize the table with NaN values
TCM = array2table(NaN(numTaxa + 1, numPermutations));

% Create row names
rowNames = strcat('coefficient_', arrayfun(@num2str, 1:numTaxa, 'UniformOutput', false));
rowNames{end+1} = 'r^2';

% Assign row names to the table
TCM.Properties.RowNames = rowNames;

% Create column names
colNames = strcat('Assemblage_', arrayfun(@num2str, 1:numPermutations, 'UniformOutput', false));

% Assign column names to the table
TCM.Properties.VariableNames = colNames;
end

function [trainingMethod, varargin, addTestData] = determineMethod(regressionMethod, settings)
if strcmp(regressionMethod, 'EQO')
    trainingMethod = @runEQO;
    varargin = {};
    addTestData = false;
else
    varargin = {regressionMethod, settings};
    trainingMethod = @runRegression;
    if ~strcmp(regressionMethod, 'OLS')
        addTestData = true; % Set the flag to true if testData and testOutput need to be added
    else
        addTestData = false;
    end
end
end

function importanceValues = getImportanceValues(TCM)
% Calculate the cumulative R^2 for each taxon
importanceValues = calculateCumulativeR2(TCM);
end

function normalizedCumulativeR2 = calculateCumulativeR2(TCM)
numTaxa = height(TCM) - 1;
numAssemblage = width(TCM);
cumulativeR2 = zeros(numTaxa, 1);
for taxonIdx = 1:numTaxa
    for permIdx = 1:numAssemblage
        cumulativeR2(taxonIdx) = cumulativeR2(taxonIdx) + TCM{taxonIdx, permIdx} * TCM{end, permIdx};
    end
end

% Normalize to the maximal importance
normalizedCumulativeR2 = cumulativeR2./max(cumulativeR2);
end