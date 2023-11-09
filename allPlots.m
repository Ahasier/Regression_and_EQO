regressionMethod = 'EQO';

numberOfTaxaInAGroup_list = 5:5:50;
numSamples_list = 10:10:200;

% Initialize all neccessary parametters from configurations files
[numPermutations, phylogenyDependency, noiseLevel, meshGrid, settings, fullIdentifier] = initializations(regressionMethod);

% Plot heatmaps
plotHeatmapsForRegressionResults(settings, regressionMethod, fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid)

%% Helper functions
function plotHeatmapsForRegressionResults(settings, regressionMethod, fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid)
fullIdentifier2 = [fullIdentifier, '_aicRegularized'];
EQOMap = loadMapFromFile('EQO', fullIdentifier2, numberOfTaxaInAGroup_list, numSamples_list, meshGrid);

% plot heatmaps
figure();
h = plotMap(numberOfTaxaInAGroup_list, numSamples_list, EQOMap, 'EQO');
saveas(gcf,['results/plots/Heatmap',fullIdentifier2,'.jpg'])
end

function Map = loadMapFromFile(regressionMethod, fullIdentifier, numberOfTaxaInAGroup_list, numSamples_list, meshGrid)
Map = accessToAccMap('results/Accuracy/', regressionMethod, fullIdentifier, 'Load');
Map = Map(numberOfTaxaInAGroup_list/meshGrid.TaxaGroup, numSamples_list/meshGrid.Samples);
end

function resultMap = accessToAccMap(resultsPath, regressionMethod, fullIdentifier, saveOrLoad)
filename = [resultsPath, 'Acc', regressionMethod, fullIdentifier, '.csv'];
if strcmp(saveOrLoad, 'Load')
    resultMap = csvread(filename);
end
end