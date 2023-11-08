function saveAccuracyToFile(regressionMethod)
% Initialize paths
initializePaths();

% Set parameters from JSON file
paramsFilename = 'configurations/basicParams.json';
[~, ~, ~, meshGrid] = setParams(paramsFilename);

% Define global variable paths for where to load data or store results
global paths
paths = SetPathsForDataAndResults('data', 'results', 'betaResults','accuracyResults', 'tcmResults');

% Set other parameters using setOptionsAndNames function
[settings, fullIdentifier] = setOptionsAndNames();

numberOfTaxaInAGroup_list = 5:5:50;
numSamples_list = 10:10:200;

% Construct the filename where accuracy results will be saved
accuracyFilename = [paths.accuracyResultsPath, 'Acc', regressionMethod, fullIdentifier, '.csv'];

if ~usingRealData(settings)
    for n = 1:length(numberOfTaxaInAGroup_list)
        numberOfTaxaInAGroup = numberOfTaxaInAGroup_list(n);
        for m = 1:length(numSamples_list)
            numSamples = numSamples_list(m);
            betaResultsFilePath = [paths.betaResultsPath, 'Betas_', regressionMethod, fullIdentifier, '_K', num2str(numberOfTaxaInAGroup), '_nSpl', num2str(numSamples), '.mat'];
            if isfile(betaResultsFilePath)
                load(betaResultsFilePath, 'results');
                
                % Save accuracy data to .csv file
                % Load or initialize existing accuracy results data
                existingData = loadOrInitializeAccuracyResultsFile(accuracyFilename);
                
                % Compute indices to determine where the accuracy result will be stored
                [index1, index2] = indicesOfAccuracyMatrixElement(numberOfTaxaInAGroup, numSamples, meshGrid);
                
                % Check if computed indices are positive integers
                if index1 <= 0 || floor(index1) ~= index1 || index2 <= 0 || floor(index2) ~= index2
                    error('Invalid indices computed: indices must be positive integers.');
                    % Check if the location specified by the indices is empty or NaN
                else
                    % Update the specified location with the new accuracy value
                    existingData(index1, index2) = results.accuracy;
                    % Save the updated dataset back to the CSV file
                    csvwrite(accuracyFilename, existingData);
                end
            end
        end
    end
end
end