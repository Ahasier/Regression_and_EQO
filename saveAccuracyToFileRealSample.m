function saveAccuracyToFileRealSample(regressionMethod, betaEps, realAbd, usePhylogeny)
% Initialize all neccessary parametters from configurations files
[~, ~, ~, meshGrid, settings] = initializations(regressionMethod);

settings.BetaEps = betaEps;

if ~strcmp(realAbd, 'EMPTY_ARRAY')
    settings.RealAbd = realAbd;
    settings.index = [];
end

if ~strcmp(usePhylogeny, 'EMPTY_ARRAY')
    settings.usePhylogeny = usePhylogeny;
end

% Create a full identifier string based on the options
fullIdentifier1 = createIdentifier(settings);

index_list = 1:25;
numSamples_list = [10 20 30 40 50 60 70 80 90 100 110 120 130];

global paths

if ~usingRealData(settings)
    for n = 1:length(index_list)
        index = index_list(n);
        settings.index = index;
        
        % Create a full identifier string based on the options
        fullIdentifier2 = createIdentifier(settings);
        
        % Construct the filename where accuracy results will be saved
        accuracyFilename = [paths.accuracyResultsPath, 'Acc', regressionMethod, fullIdentifier1, '.csv'];
        
        for m = 1:length(numSamples_list)
            numSamples = numSamples_list(m);
            betaResultsFilePath = [paths.betaResultsPath, '/Betas_', regressionMethod, fullIdentifier2, '_nSpl', num2str(numSamples), '.mat'];
            if isfile(betaResultsFilePath)
                load(betaResultsFilePath, 'results');
                
                % Save accuracy data to .csv file
                % Load or initialize existing accuracy results data
                existingData = loadOrInitializeAccuracyResultsFile(accuracyFilename);
                
                % Compute indices to determine where the accuracy result will be stored
                [index1, index2] = indicesOfAccuracyMatrix(index, numSamples, meshGrid);
                
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

function [index1, index2] = indicesOfAccuracyMatrix(index, numSamples, meshGrid)
index1 = index;
index2 = numSamples / meshGrid.Samples;
if numSamples == 136
    index2 = 14;
end
end