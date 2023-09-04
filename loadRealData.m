function [taxaAbundance, observableOutput] = loadRealData(numberOfSamples)
% LOADREALDATA loads real microbial abundance data and extracts relevant observable
% outputs (e.g., nitrate levels) from the metadata.
%
% INPUT:
%   numberOfSamples: Number of samples for which data needs to be loaded.
%
% OUTPUTS:
%   taxaAbundance: Matrix containing taxa abundance data (samples x taxa).
%   observableOutput: Vector containing observed functional outputs (samples x 1).

% Paths and constants
[pathToData, ~, ~, pathToAbundanceAndFunctionalOutputs] = SetPathsForDataAndResults();

% Load tree data and the metadata data
treeData = loadRealTreeData(pathToData);
metaTaraData = loadMetaTaraData(pathToData);

% Extract Nitrate (NO3) data and initial taxa abundance
nitrateData = extractNitrateData(metaTaraData);
initialTaxaAbundance = treeData.N1;

% Create subsets of data and save them
taxaAbundance = createAndSaveSubsetsOfData(nitrateData, initialTaxaAbundance, numberOfSamples, pathToAbundanceAndFunctionalOutputs);
observableOutput = nitrateData;
end

%% Helper functions
% Load tree data from a file
function treeData = loadRealTreeData(pathToData)
filePath = [pathToData, 'tree100taxaReal.mat'];
treeData = load(filePath);
end

% Extract nitrate data column (NO3) from Meta Tara data table
function nitrateData = extractNitrateData(metaTaraData)
nitrateData = metaTaraData{:,'NO3'};
end

% Create a random subset of data based on the number of samples in the subset,
% save the subset to files, and return this subset of taxa abundance data
function taxaAbundance = createAndSaveSubsetsOfData(nitrateData, initialTaxaAbundance, numSamples, pathToNAndX)
sampleIndices = randperm(length(nitrateData), numSamples);

ObservableOutput = nitrateData(sampleIndices);
taxaAbundance = initialTaxaAbundance(:, sampleIndices)';

saveDataToFile(taxaAbundance, [pathToNAndX, 'RealAbd_tree100_nSpl', num2str(numSamples), '.csv']);
saveDataToFile(ObservableOutput, [pathToNAndX, 'RealX_tree100_nSpl', num2str(numSamples), '.csv']);
end

% Save given data to a CSV file at the specified file path
function saveDataToFile(data, filePath)
writematrix(data, filePath);
end