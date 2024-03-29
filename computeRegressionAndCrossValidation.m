function [TCM, importanceValues, crossValidatedCoefficients] = computeRegressionAndCrossValidation(abundanceData, functionalOutput, numPermutations, regressionMethod, settings, varargin)
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
TCM = getTCM(abundanceData, functionalOutput, numPermutations, regressionMethod, settings, varargin{:});

% Check if threshold method is 'cv'
if ischar(settings.Threshold) && strcmp(settings.Threshold, 'cv')
    % If already used optimal threshold through cross-validation, no longer
    % need to go through extra crpss-validation step.
    importanceValues = [];
    crossValidatedCoefficients = median(TCM{1:end-1, :}, 2);
else
    % 2. Get importance values from TCM
    importanceValues = getImportanceValues(TCM);
    
    % 3. Get the cross-validated coefficients
    crossValidatedCoefficients = pickImportanceValue(importanceValues);
%     crossValidatedCoefficients = selectTopKTaxa(importanceValues, optimalGroupSize);
end
end

%% Helper functions
function TCM = getTCM(abundanceData, functionalOutput, numPermutations, regressionMethod, settings, varargin)
% Initialize TCM and other variables
numTaxa = size(abundanceData, 2);
TCM = initializeTCM(numTaxa, numPermutations);
aicValues = zeros(1, numPermutations);
groupSizes = zeros(1, numPermutations);

% Determine which method to use based on settings.Threshold
[trainingMethod, varargToMethod, addTestData] = determineMethod(regressionMethod, settings, varargin{:});
lenVararg = length(varargToMethod);

for i = 1:numPermutations
    % Step (a): Split data
    [trainingData, testData, trainingOutput, testOutput] = splitData(abundanceData, functionalOutput);
    
    % If incorporating phylogeny, using the grouped training abundance data
    if useExtraFeatures(settings)
        extraPhyloVars = varargin{end};
        trainingData = groupAbundanceData(trainingData, extraPhyloVars.numBranches, extraPhyloVars.addedLeaves, extraPhyloVars.Idx);
        if addTestData
            testData = groupAbundanceData(testData, extraPhyloVars.numBranches, extraPhyloVars.addedLeaves, extraPhyloVars.Idx);
        end
    end
    
    % Check if testData and testOutput need to be added to varargin
    if addTestData
        [varargToMethod{lenVararg + 1:lenVararg + 2}] = deal(testData, testOutput);
    end
    
    % Step (b): Depending on the threshold value, use regressions or EQO to train on training data
    [coefficients, aicValues(i), groupSizes(i)] = trainingMethod(trainingData, trainingOutput, varargToMethod{:});
    
    % Step (c): Calculate the out of sample R^2 on test subset
    if useExtraFeatures(settings) && addTestData
        testData = testData(:, 1:extraPhyloVars.numTaxa);
    end
    R2OutSamples = computeRSquared(testData, testOutput, coefficients);
    % referenceR2s(i) = computeRSquared(testData, testOutput, zeros(size(coefficients, 1), size(coefficients, 2)));
    
    % Store results
    TCM{:, i} = [coefficients; R2OutSamples];
end

% if sum(TCM{end, :} > 0) == 0
%     TCM{end, :} = TCM{end, :} - referenceR2s;
% end

% Determine the optimal group size
% optimalGroupSize = minimizeAIC(aicValues, groupSizes, numTaxa);
% if sum(TCM{end, :} > 0) > 0
%     optimalGroupSize = median(groupSizes(TCM{end, :} > 0));
% else
%     optimalGroupSize = median(groupSizes);
% end
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

function [trainingMethod, varargToMethod, addTestData] = determineMethod(regressionMethod, settings, varargin)
if strcmp(regressionMethod, 'EQO')
    numberOfTaxaInAGroup = varargin{1};
    trainingMethod = @runEQO;
    varargToMethod = {numberOfTaxaInAGroup, settings};
    addTestData = false;
else
    varargToMethod = {regressionMethod, settings};
    trainingMethod = @runRegression;
    
    % If doing regressions other than OLS, add testData and testOutput for
    % the cross-validation of the lambda parameter
    if ~strcmp(regressionMethod, 'OLS')
        addTestData = true;
    else
        addTestData = false;
    end
    
    % Check if phylogenetic information needs to be incorporated. If so,
    % pass extra variables related to that into trainingMethod
    if isfield(settings, 'usePhylogeny') && strcmp(settings.usePhylogeny, 'On')
        extraPhyloVars = varargin{1};
        [varargToMethod{end + 1}] = extraPhyloVars;
    end
end
end

function optimalGroupSize = minimizeAIC(aicValues, groupSizes, numTaxa)
allAicValues = zeros(numTaxa, 1);
for n = 1:numTaxa
    idxGroupSizeN = (groupSizes == n);
    allAicValues(n) = mean(aicValues(idxGroupSizeN));
end
[~, optimalGroupSize] = min(allAicValues(1:end - 10));
end

function importanceValues = getImportanceValues(TCM)
numTaxa = height(TCM) - 1;
numAssemblage = width(TCM);
cumulativeR2 = zeros(numTaxa, 1);

r2 = TCM{end, :};
r2(r2 < 0) = 0;

% Calculate the cumulative R^2 for each taxon
for taxonIdx = 1:numTaxa
    for permIdx = 1:numAssemblage
        cumulativeR2(taxonIdx) = cumulativeR2(taxonIdx) + TCM{taxonIdx, permIdx} * r2(permIdx);
    end
end

% Normalize to the maximal importance
importanceValues = cumulativeR2./max(cumulativeR2);
% importanceValues = (cumulativeR2 - min(cumulativeR2))./(max(cumulativeR2) - min(cumulativeR2));
end