function saveResults(results, resultsPath, betaResultsPath, regressionMethod, fullIdentifier, numberOfTaxaInAGroup, numSamples, accuracy, meshGrid)
% SAVERESULTS Save regression results to CSV and MAT files.
%
% INPUTS:
%   abundanceData: Matrix of abundance data.
%   functionalOuput: Values of functional outputs from different samples.
%   syntheticCoefficients: Synthetic coefficients generated.
%   resultsPath: Path to save the .CSV accuracy results.
%   betaResultsPath: Path to save the .MAT results.
%   regressionMethod: Regression method used.
%   fullIdentifier: A string indicates indentific information of the
%   settings.
%   numberOfTaxaInAGroup: Number of taxa selected in a functional group
%   when generating synthetic data.
%   numSamples: Number of total samples in the regression model.
%   estimatedCoefficients: Estimated coefficients recovered by regressions.
%   crossValidatedCoefficients: Optimal estimated coefficients selected
%   through cross-validation.
%   squaredError: Squared error of the regression model out of samples.
%   R2OutSample: R^2 correlation coefficients out of samples.
%   R2inSample: R^2 correlation coefficients in samples.
%   accuracy: Estimation accuracy of the regression model.
% 
% OUTPUTS:
%   No direct output; instead, it writes CSV and MAT files to disk.

% Construct the filename where accuracy results will be saved
filename = [resultsPath, 'Acc', regressionMethod, fullIdentifier, '.csv'];

if usingRealData(settings)
    % If using real data, skip this step
else
    % Check if data needs update
    if dataNeedsUpdate(filename, numberOfTaxaInAGroup, numSamples)
        % Compute indices to determine where the accuracy result will be stored
        [index1, index2] = indicesOfAccuracyMatrixElement(numberOfTaxaInAGroup, numSamples, meshGrid);
        
        % Load or initialize existing accuracy results data
        existingData = loadOrInitializeAccuracyResultsFile(filename);
        
        % Check if computed indices are positive integers
        if index1 <= 0 || floor(index1) ~= index1 || index2 <= 0 || floor(index2) ~= index2
            error('Invalid indices computed: indices must be positive integers.');
            % Check if the location specified by the indices is empty or NaN
        else
            % Update the specified location with the new accuracy value
            existingData(index1, index2) = accuracy;
            % Save the updated dataset back to the CSV file
            csvwrite(filename, existingData);
        end
    end
end

% Construct the file path where other simulation results will be saved
saveFilePath = [betaResultsPath, 'Betas_', regressionMethod, fullIdentifier, '_K', num2str(numberOfTaxaInAGroup), '_nSpl', num2str(numSamples), '.mat'];

% Save various results to a MAT file for later analysis
save(saveFilePath, results);
end

function doesDataNeedUpdate = dataNeedsUpdate(filename, numberOfTaxaInAGroup, numSamples)
% Load or initialize existing accuracy results data
existingData = loadOrInitializeAccuracyResultsFile(filename);

% Compute indices to determine where the accuracy result will be stored
[index1, index2] = indicesOfAccuracyMatrixElement(numberOfTaxaInAGroup, numSamples);

% Determine if data needs update
if existingData(index1, index2) == 0 || isnan(existingData(index1, index2))
    doesDataNeedUpdate = 1;
else
    doesDataNeedUpdate = 0;
end
end