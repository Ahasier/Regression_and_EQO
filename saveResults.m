function saveResults(results, regressionMethod, fullIdentifier, numberOfTaxaInAGroup, numSamples, meshGrid)
% SAVERESULTS Save regression results to CSV and MAT files.
%
% INPUTS:
%   results: The results of the regression analysis.
%   resultsPath: Path to save the .CSV accuracy results.
%   betaResultsPath: Path to save the .MAT results related to regression coefficients.
%   regressionMethod: Regression method used for the analysis.
%   fullIdentifier: A string indicating specific identification information of the settings.
%   numberOfTaxaInAGroup: Number of taxa selected in a functional group (relevant for synthetic data).
%   numSamples: Number of total samples used in the regression model.
%   meshGrid: Grid used for visualization or further analysis (if applicable).
%
% OUTPUTS:
%   Results are saved to files; no direct outputs from the function.

% Recall paths for results storage
global paths

% Construct the filename where accuracy results will be saved
filename = [paths.resultsPath, 'Acc', regressionMethod, fullIdentifier, '.csv'];

if usingRealData(settings)
    % If using real data, skip this step
else
    % Load or initialize existing accuracy results data
    existingData = loadOrInitializeAccuracyResultsFile(filename);
    
    % Compute indices to determine where the accuracy result will be stored
    [index1, index2] = indicesOfAccuracyMatrixElement(numberOfTaxaInAGroup, numSamples, meshGrid);
    
    % Check if data needs update
    if dataNeedsUpdate(filename, index1, index2)
        % Check if computed indices are positive integers
        if index1 <= 0 || floor(index1) ~= index1 || index2 <= 0 || floor(index2) ~= index2
            error('Invalid indices computed: indices must be positive integers.');
            % Check if the location specified by the indices is empty or NaN
        else
            % Update the specified location with the new accuracy value
            existingData(index1, index2) = results.accuracy;
            % Save the updated dataset back to the CSV file
            csvwrite(filename, existingData);
        end
    end
end

% Construct the file path where other simulation results will be saved
saveFilePath = [paths.betaResultsPath, 'Betas_', regressionMethod, fullIdentifier, '_K', num2str(numberOfTaxaInAGroup), '_nSpl', num2str(numSamples), '.mat'];

% Save various results to a MAT file for later analysis
save(saveFilePath, 'results');
end

%% Helper functions
function [index1, index2] = indicesOfAccuracyMatrixElement(numberOfTaxaInAGroup, numSamples, meshGrid)
index1 = numberOfTaxaInAGroup / meshGrid.TaxaGroup;
index2 = numSamples / meshGrid.Samples;
end

function doesDataNeedUpdate = dataNeedsUpdate(existingData, index1, index2)
% Determine if data needs update
if size(existingData,1) >= index1 && size(existingData,2) >= index2
    if existingData(index1, index2) == 0 || isnan(existingData(index1, index2))
        doesDataNeedUpdate = 1;
    else
        doesDataNeedUpdate = 0;
    end
else
    doesDataNeedUpdate = 1;
end
end